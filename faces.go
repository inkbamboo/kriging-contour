package kriging_contour

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"strings"
	"sync"

	"github.com/duke-git/lancet/v2/maputil"
	"github.com/panjf2000/ants/v2"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	"github.com/paulmach/orb/planar"
)

// ============================================================
//  等值面数据结构
// ============================================================

// contourEntry 等值线切割条目，用于等值面递归分割算法。
// 包含等值线的几何信息、级别以及用于方向判断的元数据。
type contourEntry struct {
	geom     orb.Geometry // 几何体：orb.Polygon（闭合环）或 orb.LineString（开线）
	level    float64      // 等值线级别值
	isClosed bool         // 是否为闭合等值线
	highSide string       // 高值侧方向标识：closed 时为 "inside"/"outside"，open 时为 "left"/"right"
	coords   []orb.Point  // 原始坐标点，用于方向判断和共边检测
}

// intersectionPt 等值线与多边形边界环的交点信息。
type intersectionPt struct {
	pt      orb.Point // 交点坐标
	ringIdx int       // 交点所在环边的段索引
	lineIdx int       // 交点在等值线上的段索引
	dist    float64   // 沿等值线的累计距离
}

// ============================================================
//  GenerateFaces 等值面生成（三步递归分割算法）
// ============================================================

// GenerateFaces 基于等值线 features 生成等值面（Polygon）features。
//
// 流程：
//  1. 解析等值线：从 features 中提取闭合环和开线，按面积排序
//  2. 遍历分割：用每条等值线依次切割边界多边形及所有子多边形
//  3. 去重合并：形状去重后并发计算每个面的 level 值
//
// 参数:
//   - features: 等值线 GeoJSON Feature 列表（来自 GenerateLines），
//     每个 feature 的 properties 含 level / is_closed / high_value_side
//   - boundary: 区域边界多边形
//   - opt: 等值线参数
//
// 返回: 等值面 GeoJSON Feature 列表，properties 含 type / level 字段
func GenerateFaces(features []*geojson.Feature, boundary orb.Polygon, opt ContourOption) []*geojson.Feature {
	if len(features) == 0 {
		return nil
	}

	boundaryPoly := closeBoundaryPolygon(boundary)

	// 1. 解析等值线为切割条目
	entries := parseContourEntries(features)
	if len(entries) == 0 {
		return nil
	}

	// 2. 分割边界多边形
	rawPolys := splitPolys(boundaryPoly, entries)
	fmt.Printf("  分割后碎片数: %d\n", len(rawPolys))

	// 3. 多边形去重（仅比较形状）
	dedupRaw := deduplicatePolygons(rawPolys)

	fmt.Printf("  去重后面数: %d\n", len(dedupRaw))

	// 4. 转换为 polyWithLevel，level 稍后计算
	mergedPolys := make([]polyWithLevel, len(dedupRaw))
	for i, poly := range dedupRaw {
		mergedPolys[i] = polyWithLevel{poly: poly}
	}
	fmt.Printf("  合并后面数: %d\n", len(mergedPolys))

	// 5. 为每个多边形计算 level
	var wg sync.WaitGroup
	p, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(body interface{}) {
		defer wg.Done()
		idx := body.(int)
		mergedPolys[idx].level = computePolygonLevel(mergedPolys[idx].poly, entries, opt)

	})
	defer p.Release()
	// 逐层提取等值线：为每个 level 单独渲染到记录 canvas
	for idx := range mergedPolys {
		wg.Add(1)
		_ = p.Invoke(idx)
	}
	wg.Wait()
	// 6. 转换为 GeoJSON Feature
	var result []*geojson.Feature
	totalArea := 0.0
	for _, mp := range mergedPolys {
		polyArea := planar.Area(mp.poly)
		totalArea += polyArea
		feat := polygonToFeature(mp.poly, mp.level, opt)
		if feat != nil {
			result = append(result, feat)
		}
	}
	return result
}

// ============================================================
//  等值线解析
// ============================================================

// parseContourEntries 从 GeoJSON features 中提取等值线切割条目。
// 按几何类型分类为闭合环和开线，闭合环按面积从大到小排序，
// 返回顺序为：闭合环（大→小）+ 开线。
func parseContourEntries(features []*geojson.Feature) []contourEntry {
	var closedRings []contourEntry
	var openLines []contourEntry

	for _, feat := range features {
		props := feat.Properties
		level, _ := props["level"].(float64)
		isClosed, _ := props["is_closed"].(bool)
		highSide, _ := props["high_value_side"].(string)

		coords := extractLineCoords(feat.Geometry)
		if len(coords) < 2 {
			continue
		}

		if isClosed {
			ring := make(orb.Ring, len(coords))
			copy(ring, coords)
			if !eqPoint(ring[0], ring[len(ring)-1], 1e-9) {
				ring = append(ring, ring[0])
			}
			poly := orb.Polygon{ring}
			if planar.Area(poly) > 0.5 {
				closedRings = append(closedRings, contourEntry{
					geom:     poly,
					level:    level,
					isClosed: true,
					highSide: highSide,
					coords:   coords,
				})
			}
		} else {
			openLines = append(openLines, contourEntry{
				geom:     orb.LineString(coords),
				level:    level,
				isClosed: false,
				highSide: highSide,
				coords:   coords,
			})
		}
	}

	// 闭合环按面积从大到小排序
	sort.Slice(closedRings, func(i, j int) bool {
		return planar.Area(closedRings[i].geom.(orb.Polygon)) > planar.Area(closedRings[j].geom.(orb.Polygon))
	})

	// 切割条目顺序：闭合环（大→小） + 开线
	entries := append(closedRings, openLines...)
	fmt.Printf("  切割条目: 闭环 %d 个, 开线 %d 条\n", len(closedRings), len(openLines))
	return entries
}

// extractLineCoords 从 orb.Geometry 中提取 LineString 坐标点列表。
// 若 geometry 为 nil 或不是 LineString 类型，返回 nil。
func extractLineCoords(geom orb.Geometry) []orb.Point {
	if geom == nil {
		return nil
	}
	ls, ok := geom.(orb.LineString)
	if !ok {
		return nil
	}
	return []orb.Point(ls)
}

// splitPolys 按 entries 顺序遍历所有等值线，逐一尝试切割多边形。
// 与 recursiveSplit 不同，此版本为迭代遍历：对每条等值线依次切分所有子多边形。
func splitPolys(poly orb.Polygon, entries []contourEntry) []orb.Polygon {
	if len(entries) == 0 {
		return []orb.Polygon{poly}
	}
	if len(entries) == 1 {
		return trySplitPolygon(poly, entries[0])
	}
	polyList := trySplitPolygon(poly, entries[0])
	for _, entry := range entries[1:] {
		var tempPolyList []orb.Polygon
		for _, subPoly := range polyList {
			tempPolyList = append(tempPolyList, trySplitPolygon(subPoly, entry)...)
		}
		polyList = tempPolyList
	}
	return polyList
}

// trySplitPolygon 尝试用单条等值线切割多边形。
// 对于闭合环：检查环是否在多边形内部，若是则按环切割（产生带孔洞的外部和环内部两个面）。
// 对于开线：检查线是否贯穿多边形，若是则沿线切割为两个独立多边形。
// 若无法切割（面积过小或其他条件不满足），返回原多边形。
func trySplitPolygon(poly orb.Polygon, entry contourEntry) []orb.Polygon {
	const minArea = 0.5
	if entry.isClosed {
		ring := entry.geom.(orb.Polygon)[0]
		if planar.Area(orb.Polygon{ring}) < minArea {
			return []orb.Polygon{poly}
		}
		if isRingInsidePolygon(ring, poly) {
			return splitByClosedRing(poly, ring, minArea)
		}
	} else {
		line := entry.geom.(orb.LineString)
		if len(line) < 2 {
			return []orb.Polygon{poly}
		}
		if !isOpenLineThroughPolygon(line, poly) {
			return []orb.Polygon{poly}
		}
		result := splitByOpenLine(poly, line, minArea)
		if len(result) >= 2 {
			return result
		}
	}
	return []orb.Polygon{poly}
}

// isOpenLineThroughPolygon 判断开线是否贯穿多边形。
// 贯穿的条件：
//  1. 线的首尾端点都在多边形外环边界上（容差内）
//  2. 首尾端点不重合
//  3. 线的中间部分穿过多边形内部（至少有一个中间点严格在内部）
func isOpenLineThroughPolygon(line orb.LineString, poly orb.Polygon) bool {
	if len(line) < 2 {
		return false
	}
	exterior := poly[0]
	if len(exterior) < 4 {
		return false
	}

	startPt := line[0]
	endPt := line[len(line)-1]

	// 条件1: 首尾点必须在多边形边界上
	if minDistToRing(startPt, exterior) > 0.01 ||
		minDistToRing(endPt, exterior) > 0.01 {
		return false
	}

	// 条件2: 首尾点不能重合
	if eqPoint(startPt, endPt, 0.001) {
		return false
	}

	// 条件3: 至少有一个中间点严格在多边形内部（不在孔洞中）
	hasPointInside := false
	for _, pt := range line {
		if !planar.PolygonContains(poly, pt) {
			continue
		}
		// 再确认该点不在边界上（严格内部）
		if minDistToRing(pt, exterior) > 0.01 {
			hasPointInside = true
			break
		}
	}
	if !hasPointInside {
		return false
	}

	return true
}

// ============================================================
//  闭合环切割
// ============================================================

// isRingInsidePolygon 判断闭合环是否完全位于多边形内部。
// 检查环的所有顶点及中心点是否都在多边形内（不在边界和孔洞中）。
func isRingInsidePolygon(ring orb.Ring, poly orb.Polygon) bool {
	if len(ring) < 3 {
		return false
	}
	// 检查所有顶点都在多边形内
	for _, pt := range ring {
		if !planar.PolygonContains(poly, pt) {
			return false
		}
	}
	// 额外检查中心点，确保不在边界上
	cx, cy := 0.0, 0.0
	for _, pt := range ring {
		cx += pt[0]
		cy += pt[1]
	}
	n := float64(len(ring))
	cx /= n
	cy /= n
	return planar.PolygonContains(poly, orb.Point{cx, cy})
}

// splitByClosedRing 用闭合环切割多边形。
// 返回两个多边形：
//   - 外部多边形：原多边形 + 新增孔洞（孔洞方向为 CW）
//   - 内部多边形：环本身作为实心面（外环方向为 CCW）
//
// 遵循右手定则：外环 CCW，孔洞 CW。
func splitByClosedRing(poly orb.Polygon, ring orb.Ring, minArea float64) []orb.Polygon {
	// 外部多边形：原多边形 + 新增孔洞
	outside := make(orb.Polygon, 0, len(poly)+1)
	outside = append(outside, poly[0]) // 外环不变
	for i := 1; i < len(poly); i++ {
		outside = append(outside, poly[i]) // 保留原有孔洞
	}

	// 孔洞必须为 CW（有符号面积为负）
	holeRing := make(orb.Ring, len(ring))
	copy(holeRing, ring)
	if signedRingArea(holeRing) > 0 {
		reverseRing(holeRing)
	}
	outside = append(outside, holeRing) // 新增孔洞

	// 内部多边形：环本身作为实心多边形（外环应为 CCW）
	innerRing := make(orb.Ring, len(ring))
	copy(innerRing, ring)
	if signedRingArea([]orb.Point(innerRing)) < 0 {
		reverseRing(innerRing)
	}
	inside := orb.Polygon{innerRing}

	var result []orb.Polygon
	if planar.Area(outside) > minArea {
		result = append(result, outside)
	}
	if planar.Area(inside) > minArea {
		result = append(result, inside)
	}
	if len(result) < 2 {
		return []orb.Polygon{poly}
	}
	return result
}

// reverseRing 原地反转环的顶点顺序（双指针法）。
func reverseRing(ring orb.Ring) {
	for i, j := 0, len(ring)-1; i < j; i, j = i+1, j-1 {
		ring[i], ring[j] = ring[j], ring[i]
	}
}

// minDistToRing 计算点到环的所有边的最短距离。
func minDistToRing(pt orb.Point, ring orb.Ring) float64 {
	minD := math.MaxFloat64
	for i := 0; i < len(ring)-1; i++ {
		d := pointToSegmentDist(pt, ring[i], ring[i+1])
		if d < minD {
			minD = d
		}
	}
	return minD
}

// pointToSegmentDist 计算点到线段的最短欧氏距离。
// 通过将点投影到线段上并钳制 t 参数到 [0, 1] 实现。
func pointToSegmentDist(p, a, b orb.Point) float64 {
	dx, dy := b[0]-a[0], b[1]-a[1]
	lenSq := dx*dx + dy*dy
	if lenSq < 1e-15 {
		return math.Hypot(p[0]-a[0], p[1]-a[1])
	}
	t := ((p[0]-a[0])*dx + (p[1]-a[1])*dy) / lenSq
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	proj := orb.Point{a[0] + t*dx, a[1] + t*dy}
	return math.Hypot(p[0]-proj[0], p[1]-proj[1])
}

// lineLength 计算 LineString 的总长度（所有线段长度之和）。
func lineLength(line orb.LineString) float64 {
	total := 0.0
	for i := 0; i < len(line)-1; i++ {
		total += math.Hypot(line[i+1][0]-line[i][0], line[i+1][1]-line[i][1])
	}
	return total
}

// splitByOpenLine 用开线将多边形切割为两个独立多边形。
//
// 算法：
//  1. 验证线首尾点均在多边形边界上
//  2. 将线的首尾投影到环上作为交点
//  3. 将环在两个交点处拆分，与线的正向/反向段拼接成两个闭合多边形
//  4. 将原多边形的孔洞分配给包含它的新多边形
func splitByOpenLine(poly orb.Polygon, line orb.LineString, minArea float64) []orb.Polygon {
	exterior := poly[0]

	// 验证线首尾点必须在多边形边界上
	if minDistToRing(line[0], exterior) > 0.01 || minDistToRing(line[len(line)-1], exterior) > 0.01 {
		return []orb.Polygon{poly}
	}

	// 线的首尾点就是交点
	firstPt, _ := nearestPointOnRing(line[0], poly[0])
	lastPt, _ := nearestPointOnRing(line[len(line)-1], poly[0])
	first := intersectionPt{pt: firstPt, lineIdx: 0, dist: 0}
	last := intersectionPt{pt: lastPt, lineIdx: len(line) - 2, dist: lineLength(line)}

	if eqPoint(first.pt, last.pt, 0.001) {
		return []orb.Polygon{poly}
	}

	// 沿等值线提取 first→last 的线段（整条线）
	cutForward := extractLinePortion(line, first, last)
	// 反向线段
	cutBackward := reverseLineString(cutForward)

	arcAB, arcBA := splitRingAtPoints(poly[0], first.pt, last.pt)

	poly1 := buildPolygonFromSegments(cutForward, arcBA)
	poly2 := buildPolygonFromSegments(cutBackward, arcAB)

	// 将原多边形的孔洞分配给包含它的新多边形
	for i := 1; i < len(poly); i++ {
		hole := poly[i]
		testPt := hole[0]
		if planar.PolygonContains(poly1, testPt) {
			poly1 = append(poly1, hole)
		} else if planar.PolygonContains(poly2, testPt) {
			poly2 = append(poly2, hole)
		}
	}

	var result []orb.Polygon
	if planar.Area(poly1) > minArea {
		result = append(result, poly1)
	}
	if planar.Area(poly2) > minArea {
		result = append(result, poly2)
	}

	if len(result) < 2 {
		return []orb.Polygon{poly}
	}
	return result
}

// nearestPointOnRing 计算点 p 到环 ring 的最近投影点及距离。
// 遍历环的所有边，对每条边调用 closestPointOnSegment，返回最短距离对应的投影点。
func nearestPointOnRing(pt orb.Point, ring orb.Ring) (orb.Point, float64) {
	if len(ring) < 2 {
		return pt, 0
	}

	bestPt := orb.Point{}
	minDist := math.MaxFloat64

	for i := 0; i < len(ring)-1; i++ {
		a, b := ring[i], ring[i+1]
		closest := closestPointOnSegment(pt, a, b)
		d := math.Hypot(pt[0]-closest[0], pt[1]-closest[1])
		if d < minDist {
			minDist = d
			bestPt = closest
		}
	}

	return bestPt, minDist
}

// splitRingAtPoints 使用环上两个点 A B 将环切割成两条线。
// A、B 不要求是环的节点，可以是环上的任意点（容差 0.01）。
// 返回两条 orb.LineString: A→B 沿环的弧段, B→A 沿环的弧段（绕回）。
// 若任一点不在环上或两点重合，返回 nil, nil。
//
// 实现：找到 A、B 在环的哪段边上，将其作为节点插入环，
// 然后沿环提取 A→B 和 B→A 两段。
func splitRingAtPoints(ring orb.Ring, pA, pB orb.Point) (orb.LineString, orb.LineString) {
	if len(ring) < 3 {
		return nil, nil
	}
	if minDistToRing(pA, ring) > 0.01 || minDistToRing(pB, ring) > 0.01 {
		return nil, nil
	}
	if eqPoint(pA, pB, 0.001) {
		return nil, nil
	}

	// 获取独立顶点（去掉闭合重复点）
	n := len(ring)
	if eqPoint(ring[0], ring[n-1], 1e-9) {
		n--
	}
	verts := ring[:n]

	// 将 pA、pB 插入顶点列表
	expanded := insertPointsIntoVerts(verts, pA, pB)

	// 找 pA、pB 在扩展列表中的位置
	idxA := indexOfPoint(expanded, pA)
	idxB := indexOfPoint(expanded, pB)
	if idxA < 0 || idxB < 0 {
		return nil, nil
	}

	arcAB := extractArcForward(expanded, idxA, idxB)
	arcBA := extractArcForward(expanded, idxB, idxA)
	return arcAB, arcBA
}

// insertPointsIntoVerts 将 pA、pB 插入到环的顶点列表中。
// 若点已是已有顶点（容差 0.001）则不重复插入；
// 否则插入到其所在段的对应位置，同段上的多个点按沿段方向排列。
func insertPointsIntoVerts(verts []orb.Point, pA, pB orb.Point) []orb.Point {
	n := len(verts)

	type locateResult struct {
		pt     orb.Point
		segIdx int
		t      float64 // 沿段位置 [0,1]
		onVtx  bool    // 是否已是顶点
	}

	locate := func(pt orb.Point) locateResult {
		// 先检查是否已是顶点
		for i, v := range verts {
			if eqPoint(pt, v, 0.001) {
				return locateResult{pt, i, 0, true}
			}
		}
		// 找最近段及投影参数 t
		minDist := math.MaxFloat64
		bestSeg, bestT := 0, 0.0
		for i := 0; i < n; i++ {
			a, b := verts[i], verts[(i+1)%n]
			dx, dy := b[0]-a[0], b[1]-a[1]
			lenSq := dx*dx + dy*dy
			var t float64
			if lenSq < 1e-15 {
				t = 0
			} else {
				t = ((pt[0]-a[0])*dx + (pt[1]-a[1])*dy) / lenSq
				t = math.Max(0, math.Min(1, t))
			}
			proj := orb.Point{a[0] + t*dx, a[1] + t*dy}
			d := math.Hypot(pt[0]-proj[0], pt[1]-proj[1])
			if d < minDist {
				minDist = d
				bestSeg = i
				bestT = t
			}
		}
		return locateResult{pt, bestSeg, bestT, false}
	}

	locA := locate(pA)
	locB := locate(pB)

	// 按段收集需插入的点（跳过已是顶点的）
	segInserts := make([][]locateResult, n)
	for _, loc := range []locateResult{locA, locB} {
		if loc.onVtx {
			continue
		}
		segInserts[loc.segIdx] = append(segInserts[loc.segIdx], loc)
	}

	// 同段上多个点按 t 排序（沿段方向 A→B）
	for i := range segInserts {
		if len(segInserts[i]) > 1 {
			sort.Slice(segInserts[i], func(a, b int) bool {
				return segInserts[i][a].t < segInserts[i][b].t
			})
		}
	}

	// 构建扩展列表：遍历原顶点，每段后插入对应点
	result := make([]orb.Point, 0, n+2)
	for i := 0; i < n; i++ {
		result = append(result, verts[i])
		for _, ins := range segInserts[i] {
			result = append(result, ins.pt)
		}
	}

	return result
}

// indexOfPoint 在顶点列表中查找点的索引，使用容差匹配（0.001）。
// 返回索引值，未找到返回 -1。
func indexOfPoint(verts []orb.Point, pt orb.Point) int {
	for i, v := range verts {
		if eqPoint(v, pt, 0.001) {
			return i
		}
	}
	return -1
}

// extractArcForward 从顶点列表中提取从 from 到 to 的弧段。
// 沿环前进方向，支持跨末端环绕（from > to 时自动绕回）。
func extractArcForward(verts []orb.Point, from, to int) orb.LineString {
	n := len(verts)
	if from == to {
		return orb.LineString{verts[from]}
	}
	var pts []orb.Point
	if from < to {
		pts = make([]orb.Point, to-from+1)
		copy(pts, verts[from:to+1])
	} else {
		pts = make([]orb.Point, 0, n-from+to+1)
		pts = append(pts, verts[from:]...)
		pts = append(pts, verts[:to+1]...)
	}
	return orb.LineString(pts)
}

// extractLinePortion 沿等值线提取从 first 交点到 last 交点之间的线段。
// 包含首尾交点及之间的中间顶点，跳过重复点。
func extractLinePortion(line orb.LineString, first, last intersectionPt) orb.LineString {
	var pts orb.LineString
	pts = append(pts, first.pt)

	// 添加 first 和 last 之间的中间顶点
	for li := first.lineIdx + 1; li <= last.lineIdx; li++ {
		if !eqPoint(pts[len(pts)-1], line[li], 1e-9) {
			pts = append(pts, line[li])
		}
	}

	if !eqPoint(pts[len(pts)-1], last.pt, 1e-9) {
		pts = append(pts, last.pt)
	}

	return pts
}

// reverseLineString 反转 LineString 的顶点顺序，返回新的 LineString。
func reverseLineString(ls orb.LineString) orb.LineString {
	n := len(ls)
	result := make(orb.LineString, n)
	for i := 0; i < n; i++ {
		result[i] = ls[n-1-i]
	}
	return result
}

// buildPolygonFromSegments 将两条线段拼接为闭合多边形。
// seg1 和 seg2 首尾相连形成闭合环，跳过重合点，确保闭合并去除连续重复点。
func buildPolygonFromSegments(seg1, seg2 orb.LineString) orb.Polygon {
	var ring orb.Ring
	ring = append(ring, seg1...)
	// seg2 的第一个点与 seg1 的最后一个点应该重合，跳过
	if len(seg2) > 0 && len(ring) > 0 && eqPoint(ring[len(ring)-1], seg2[0], 1e-9) {
		ring = append(ring, seg2[1:]...)
	} else {
		ring = append(ring, seg2...)
	}
	// 确保闭合
	if len(ring) > 0 && !eqPoint(ring[0], ring[len(ring)-1], 1e-9) {
		ring = append(ring, ring[0])
	}
	// 去除连续重复点
	ring = cleanRing(ring)
	if len(ring) < 4 {
		return orb.Polygon{ring}
	}
	return orb.Polygon{ring}
}

// cleanRing 去除环中连续的重复顶点，保留第一个出现的点。
// 若顶点与上一个保留点距离 < 1e-9 则跳过。
func cleanRing(ring orb.Ring) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	result := orb.Ring{ring[0]}
	for i := 1; i < len(ring); i++ {
		if !eqPoint(result[len(result)-1], ring[i], 1e-9) {
			result = append(result, ring[i])
		}
	}
	return result
}

// ============================================================
//  多边形去重
// ============================================================

// polyWithLevel 多边形与 level 值的关联结构体。
type polyWithLevel struct {
	poly  orb.Polygon
	level float64
}

// deduplicatePolygons 对多边形列表进行纯形状去重。
// 通过 polygonKey 生成每个多边形的归一化标识键，相同形状只保留第一个。
func deduplicatePolygons(polys []orb.Polygon) []orb.Polygon {
	seen := make(map[string]bool)
	var result []orb.Polygon
	for _, poly := range polys {
		key := polygonKey(poly)
		if seen[key] {
			continue
		}
		seen[key] = true
		result = append(result, poly)
	}
	return result
}

// polygonKey 为多边形生成唯一标识键。
// 先将每个环归一化（旋转到最小顶点为起点），再序列化为字符串：
// 格式为 "x1,y1;x2,y2;...|x1,y1;x2,y2;..."。
func polygonKey(poly orb.Polygon) string {
	var parts []string
	for _, ring := range poly {
		normRing := normalizeRing(ring)
		var pts []string
		for _, pt := range normRing {
			pts = append(pts, fmt.Sprintf("%.6f,%.6f", pt[0], pt[1]))
		}
		parts = append(parts, strings.Join(pts, ";"))
	}
	return strings.Join(parts, "|")
}

// normalizeRing 归一化环：旋转到以最小坐标点（lexicographically）为起点，确保首尾闭合。
func normalizeRing(ring orb.Ring) orb.Ring {
	if len(ring) < 3 {
		return ring
	}
	n := len(ring) - 1 // 排除闭合重复点
	if n < 2 {
		return ring
	}

	// 找到 lexicographically 最小的点
	minIdx := 0
	for i := 1; i < n; i++ {
		if ring[i][0] < ring[minIdx][0] ||
			(math.Abs(ring[i][0]-ring[minIdx][0]) < 1e-9 && ring[i][1] < ring[minIdx][1]) {
			minIdx = i
		}
	}

	result := make(orb.Ring, 0, n+1)
	for i := 0; i < n; i++ {
		result = append(result, ring[(minIdx+i)%n])
	}
	result = append(result, result[0])
	return result
}

// calcSharedBoundaryLen 计算两个环之间共享边界的总长度。
// 通过顶点匹配检测共享边段（同向或反向），累加匹配的边段长度。
func calcSharedBoundaryLen(r1, r2 orb.Ring) float64 {
	n1, n2 := len(r1)-1, len(r2)-1
	total := 0.0

	for i := 0; i < n1; i++ {
		for j := 0; j < n2; j++ {
			// 检查 r1[i]→r1[i+1] 与 r2[j]→r2[j+1] 是否共享
			if (eqPoint(r1[i], r2[j], 0.001) && eqPoint(r1[i+1], r2[j+1], 0.001)) ||
				(eqPoint(r1[i], r2[j+1], 0.001) && eqPoint(r1[i+1], r2[j], 0.001)) {
				segLen := math.Hypot(r1[i+1][0]-r1[i][0], r1[i+1][1]-r1[i][1])
				total += segLen
			}
		}
	}
	return total
}

// ============================================================
//  多边形与线的位置关系
// ============================================================

// polygonSideOfLine 计算多边形相对于等值线的位置关系。
//
// 若多边形与等值线共边，则判断多边形在等值线的哪一侧。
//   - entry.isClosed: 返回 "inside" 或 "outside"（相对环的位置）
//   - entry 为开线: 返回 "left" 或 "right"（相对线走向的位置）
//
// hasSharedEdge 表示是否存在共边关系。
func polygonSideOfLine(poly orb.Polygon, entry contourEntry) (side string, hasSharedEdge bool) {
	if len(entry.coords) < 2 {
		return "", false
	}
	// 判断是否与外环共边
	if !sharesEdgeWithExterior(poly[0], entry.geom) {
		return "", false
	}
	repPt := representativePoint(poly)
	if entry.isClosed {
		ring := entry.geom.(orb.Polygon)[0]
		ringPoly := orb.Polygon{ring}
		if planar.PolygonContains(ringPoly, repPt) {
			return "inside", true
		}
		return "outside", true
	}

	// 对于开线，找到共享边段，根据共享边的方向判断多边形在哪一侧
	start, end, reversed, found := findSharedEdgeSegment(poly[0], entry.geom)

	if found {
		// 计算共享边的向量和长度
		dx := end[0] - start[0]
		dy := end[1] - start[1]
		len := math.Hypot(dx, dy)
		if len < 1e-12 {
			// 退化线段，无法判断方向，默认返回右侧
			return "right", true
		}
		// 取共享边上的中点
		mid := orb.Point{(start[0] + end[0]) / 2, (start[1] + end[1]) / 2}
		// 左侧法向量为 (-dy/len, dx/len)
		epsilon := 0.001
		leftPt := orb.Point{mid[0] - dy/len*epsilon, mid[1] + dx/len*epsilon}
		// 检查左侧点是否在多边形内部
		isLeft := planar.PolygonContains(poly, leftPt)
		// 如果等值线方向与环方向相反，结果翻转
		if reversed {
			isLeft = !isLeft
		}
		if isLeft {
			return "left", true
		}
		return "right", true
	}

	// 回退：使用整体走向（原始逻辑）
	return pointSideOfLine(repPt, entry.coords), true
}

// ============================================================
//
//	Level 计算
//
// ============================================================
// computePolygonLevel 基于等值线与多边形的空间关系推断多边形的 level 值。
// 遍历所有等值线，对每条与多边形共边的等值线：
//   - 若 polygonSide 与 entry.highSide 一致，level = entry.level
//   - 否则 level = entry.level - opt.ContourInterval
//
// 返回第一个匹配结果即退出。
func computePolygonLevel(poly orb.Polygon, entries []contourEntry, opt ContourOption) (level float64) {
	for _, entry := range entries {
		side, hasSharedEdge := polygonSideOfLine(poly, entry)
		if !hasSharedEdge {
			continue
		}
		if entry.highSide == side {
			return entry.level
		}
		return math.Round((entry.level-opt.ContourInterval)*1e2) / 1e2
	}
	return
}

// sharesEdgeWithExterior 判断等值线几何体是否与多边形的某个环有非退化共边。
// 使用端点匹配法：检查是否有至少一个环的线段与等值线线段共享两个端点。
func sharesEdgeWithExterior(ring orb.Ring, geom orb.Geometry) bool {
	var contourLS orb.LineString
	switch g := geom.(type) {
	case orb.Polygon:
		if len(g) == 0 {
			return false
		}
		contourLS = orb.LineString(g[0])
	case orb.LineString:
		contourLS = g
	default:
		return false
	}

	n1, n2 := len(ring)-1, len(contourLS)
	if n2 < 2 {
		return false
	}

	totalMatch := 0.0
	for i := 0; i < n1; i++ {
		for j := 0; j < n2-1; j++ {
			// 检查环段 (ring[i],ring[i+1]) 与等值线段 (contourLS[j],contourLS[j+1]) 是否端点匹配
			// 支持同向和反向
			if (eqPoint(ring[i], contourLS[j], 0.001) && eqPoint(ring[i+1], contourLS[j+1], 0.001)) ||
				(eqPoint(ring[i], contourLS[j+1], 0.001) && eqPoint(ring[i+1], contourLS[j], 0.001)) {
				segLen := math.Hypot(ring[i+1][0]-ring[i][0], ring[i+1][1]-ring[i][1])
				totalMatch += segLen
			}
		}
	}
	return totalMatch > 0.01
}

// findSharedEdgeSegment 找到多边形环与等值线几何体共享的边段。
//
// 返回值：
//   - start, end: 共享边的起点和终点（按多边形环的方向）
//   - reversed: 等值线方向与环方向是否相反
//   - found: 是否找到共享边段
func findSharedEdgeSegment(ring orb.Ring, geom orb.Geometry) (start, end orb.Point, reversed bool, found bool) {
	var contourLS orb.LineString
	switch g := geom.(type) {
	case orb.Polygon:
		if len(g) == 0 {
			return
		}
		contourLS = orb.LineString(g[0])
	case orb.LineString:
		contourLS = g
	default:
		return
	}

	n1, n2 := len(ring)-1, len(contourLS)
	if n2 < 2 {
		return
	}

	for i := 0; i < n1; i++ {
		for j := 0; j < n2-1; j++ {
			// 检查环段 (ring[i],ring[i+1]) 与等值线段 (contourLS[j],contourLS[j+1]) 是否端点匹配
			if eqPoint(ring[i], contourLS[j], 0.001) && eqPoint(ring[i+1], contourLS[j+1], 0.001) {
				// 同向匹配
				return ring[i], ring[i+1], false, true
			} else if eqPoint(ring[i], contourLS[j+1], 0.001) && eqPoint(ring[i+1], contourLS[j], 0.001) {
				// 反向匹配：等值线方向与环方向相反
				return ring[i], ring[i+1], true, true
			}
		}
	}
	return
}

// representativePoint 取多边形内部的一个保证不在边界上的代表点。
// 优先使用外环中心点，若在孔洞中则取第一条边的中垂线内侧点。
func representativePoint(poly orb.Polygon) orb.Point {
	if len(poly) == 0 || len(poly[0]) < 4 {
		return orb.Point{0, 0}
	}

	// 取外环的中心点
	cx, cy := 0.0, 0.0
	n := len(poly[0]) - 1 // 排除闭合点
	for i := 0; i < n; i++ {
		cx += poly[0][i][0]
		cy += poly[0][i][1]
	}
	cx /= float64(n)
	cy /= float64(n)

	// 验证中心点在多边形内部（不在孔洞中）
	center := orb.Point{cx, cy}
	if planar.PolygonContains(poly, center) {
		return center
	}

	// 回退：尝试取第一条边的中垂线内侧点
	if n >= 3 {
		a, b := poly[0][0], poly[0][1]
		mid := orb.Point{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2}
		dx, dy := b[0]-a[0], b[1]-a[1]
		l := math.Hypot(dx, dy)
		if l > 1e-9 {
			step := 0.001
			pt := orb.Point{mid[0] - dy/l*step, mid[1] + dx/l*step}
			if planar.PolygonContains(poly, pt) {
				return pt
			}
			pt2 := orb.Point{mid[0] + dy/l*step, mid[1] - dx/l*step}
			if planar.PolygonContains(poly, pt2) {
				return pt2
			}
		}
	}

	return center
}

// pointSideOfLine 判断点在等值线总体走向的左侧还是右侧。
// 使用二维叉积：cross > 0 为左，否则为右。
func pointSideOfLine(pt orb.Point, coords []orb.Point) string {
	if len(coords) < 2 {
		return "right"
	}
	// 使用等值线的总体方向
	sx, sy := coords[0][0], coords[0][1]
	ex, ey := coords[len(coords)-1][0], coords[len(coords)-1][1]
	dx := ex - sx
	dy := ey - sy

	if math.Abs(dx) < 1e-9 && math.Abs(dy) < 1e-9 {
		return "right"
	}

	px := pt[0] - sx
	py := pt[1] - sy
	cross := dx*py - dy*px
	if cross > 0 {
		return "left"
	}
	return "right"
}

// ============================================================
//  GeoJSON 输出
// ============================================================

// polygonToFeature 将 orb.Polygon 转换为 GeoJSON Feature。
// 对多边形坐标四舍五入（6 位小数），properties 含 type / level 及透传的 opt.Extra。
func polygonToFeature(poly orb.Polygon, level float64, opt ContourOption) *geojson.Feature {
	if len(poly) == 0 {
		return nil
	}

	// 确保外环是 CCW，孔洞是 CW
	exterior := poly[0]
	if len(exterior) < 4 {
		return nil
	}

	// 构建 orb.Polygon
	resultPoly := make(orb.Polygon, 0, len(poly))
	resultPoly = append(resultPoly, roundRing(exterior))
	for i := 1; i < len(poly); i++ {
		resultPoly = append(resultPoly, roundRing(poly[i]))
	}

	feature := geojson.NewFeature(resultPoly)
	feature.Properties = maputil.Merge(map[string]interface{}{
		"type":  opt.ContourType,
		"level": math.Round(level*100) / 100,
	}, opt.Extra)

	return feature
}

// roundRing 将环的所有顶点坐标四舍五入到 6 位小数。
func roundRing(ring orb.Ring) orb.Ring {
	result := make(orb.Ring, len(ring))
	for i, pt := range ring {
		result[i] = orb.Point{math.Round(pt[0]*1e6) / 1e6, math.Round(pt[1]*1e6) / 1e6}
	}
	return result
}
