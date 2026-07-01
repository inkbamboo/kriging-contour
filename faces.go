// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了等值面生成的算法实现，包括多边形切割、区域合并、属性计算等。
package kriging_contour

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/duke-git/lancet/v2/maputil"
	"github.com/inkbamboo/kriging-contour/internal/utils"
	"github.com/panjf2000/ants/v2"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	"github.com/paulmach/orb/planar"
)

// contourEntry 表示一条等值线条目，包含其几何对象、level 值、闭合状态、高值侧和坐标点。
// 用于等值面切割流程中的中间数据结构。
type contourEntry struct {
	geom     orb.Geometry // 几何对象（Polygon 表示闭合环，LineString 表示开放线）
	level    float64      // 等值线级别值
	isClosed bool         // 是否闭合
	highSide string       // 高值侧方向："inside"/"outside" 或 "left"/"right"
	coords   []*orb.Point // 坐标点序列
}

// intersectionPt 表示等值线与边界的一个交点。
type intersectionPt struct {
	pt      orb.Point // 交点坐标
	ringIdx int       // 交点在多边形环中的段索引
	lineIdx int       // 交点在等值线中的段索引
	dist    float64   // 交点到等值线起点的累积距离
}

// GenerateFaces 从等值线 features 生成等值面（填充多边形）。
// 这是等值面生成的主入口函数。
//
// 完整流程：
//  1. 从 features 中分离边界特征和等值线特征
//  2. 解析 contourEntry 列表（闭环保序按面积降序，开线在后）
//  3. 递归用每条等值线切割边界多边形（splitPolys）
//  4. 去重（deduplicatePolygons）
//  5. 并行计算每个面的 level（computePolygonLevel）
//  6. 转换为 GeoJSON Feature 并输出
//
// 参数：
//   - features: GenerateLines 返回的等值线 GeoJSON Feature 列表（含边界 Feature）
//   - boundary: 原始边界多边形
//   - opt: 等值线配置参数
//
// 返回等值面 GeoJSON Feature 切片。
func GenerateFaces(features []*geojson.Feature, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature {
	start := time.Now()
	fmt.Printf("[INFO] 步骤开始: 等值面生成\n")
	defer func() {
		fmt.Printf("[INFO] 步骤完成: 等值面生成 (耗时: %v)\n", time.Since(start).Round(time.Millisecond))
	}()

	if len(features) == 0 {
		fmt.Printf("[WARN] 等值线 features 为空，无法生成等值面\n")
		return nil
	}

	var contourFeatures []*geojson.Feature
	for _, f := range features {
		if f.Properties["type"] == "boundary" {
			continue
		}
		contourFeatures = append(contourFeatures, f)
	}
	fmt.Printf("[INFO] 等值线特征数: %d (不含边界)\n", len(contourFeatures))

	closeBoundaryPolygon(boundary)
	boundaryPoly := boundary

	entries := parseContourEntries(contourFeatures)
	if len(entries) == 0 {
		fmt.Printf("[WARN] 解析等值线条目为空，无法生成等值面\n")
		return nil
	}

	rawPolys := splitPolys(boundaryPoly, entries)
	fmt.Printf("[INFO] 分割后碎片数: %d\n", len(rawPolys))
	if len(rawPolys) == 0 {
		fmt.Printf("[WARN] 分割后无碎片，返回空等值面\n")
		return nil
	}

	dedupRaw := deduplicatePolygons(rawPolys)
	fmt.Printf("[INFO] 去重后面数: %d\n", len(dedupRaw))
	if len(dedupRaw) == 0 {
		fmt.Printf("[WARN] 去重后无有效面，返回空等值面\n")
		return nil
	}

	mergedPolys := make([]polyWithLevel, len(dedupRaw))
	for i, poly := range dedupRaw {
		mergedPolys[i] = polyWithLevel{poly: poly}
	}
	fmt.Printf("[INFO] 合并后面数: %d\n", len(mergedPolys))

	levelStart := time.Now()
	fmt.Printf("[INFO] 步骤开始: 等值面 level 计算\n")
	var wg sync.WaitGroup
	// 并行计算每个面的 level
	p, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(body interface{}) {
		defer wg.Done()
		idx := body.(int)
		mergedPolys[idx].level = computePolygonLevel(mergedPolys[idx].poly, entries, opt)
	})
	defer p.Release()
	for idx := range mergedPolys {
		wg.Add(1)
		_ = p.Invoke(idx)
	}
	wg.Wait()
	fmt.Printf("[INFO] 步骤完成: 等值面 level 计算 (耗时: %v)\n", time.Since(levelStart).Round(time.Millisecond))

	var result []*geojson.Feature
	totalArea := 0.0
	emptyCount := 0
	for _, mp := range mergedPolys {
		if len(mp.poly) == 0 || len(mp.poly[0]) < 4 {
			emptyCount++
			continue
		}
		polyArea := planar.Area(mp.poly)
		totalArea += polyArea
		feat := polygonToFeature(mp.poly, mp.level, opt)
		if feat != nil {
			result = append(result, feat)
		}
	}
	if emptyCount > 0 {
		fmt.Printf("[DEBUG] 跳过了 %d 个空面\n", emptyCount)
	}
	fmt.Printf("[INFO] 最终等值面数: %d, 总面积: %.2f\n", len(result), totalArea)
	return result
}

// parseContourEntries 从等值线 GeoJSON Feature 列表中解析 contourEntry 列表。
// 闭环按面积降序排列（大环优先处理），开线排在闭环之后。
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
			ring := orb.Ring(toOrbPoints(coords))
			if !utils.EqPoint(ring[0], ring[len(ring)-1], 1e-9) {
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
				geom:     orb.LineString(toOrbPoints(coords)),
				level:    level,
				isClosed: false,
				highSide: highSide,
				coords:   coords,
			})
		}
	}

	sort.Slice(closedRings, func(i, j int) bool {
		return planar.Area(closedRings[i].geom.(orb.Polygon)) > planar.Area(closedRings[j].geom.(orb.Polygon))
	})

	entries := append(closedRings, openLines...)
	fmt.Printf("[DEBUG] 切割条目: 闭环 %d 个, 开线 %d 条\n", len(closedRings), len(openLines))
	return entries
}

// extractLineCoords 从 orb.Geometry 中提取坐标点序列（仅支持 LineString 类型）。
func extractLineCoords(geom orb.Geometry) []*orb.Point {
	if geom == nil {
		return nil
	}
	ls, ok := geom.(orb.LineString)
	if !ok {
		return nil
	}
	return fromOrbPoints(ls)
}

// splitPolys 递归用所有 contourEntry 切割多边形。
// 第一个 entry 切割原始多边形，后续 entry 逐步切割结果列表。
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

// trySplitPolygon 尝试用单条 contourEntry 切割一个多边形。
// 闭合环使用 splitByClosedRing（形成外环+孔洞结构），
// 开放线使用 splitByOpenLine（沿线段一分为二）。
// 若切割条件不满足或面积过小，返回原多边形。
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

// isOpenLineThroughPolygon 判断开放线是否穿过给定的多边形。
// 要求：线的起止点均在边界上，且线内部有点在多边形内部（不紧贴边界）。
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

	if minDistToRing(startPt, exterior) > 0.01 ||
		minDistToRing(endPt, exterior) > 0.01 {
		return false
	}

	if utils.EqPoint(startPt, endPt, 0.001) {
		return false
	}

	hasPointInside := false
	for _, pt := range line {
		if !planar.PolygonContains(poly, pt) {
			continue
		}
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

// isRingInsidePolygon 判断环是否完全位于多边形内部。
// 检查环的所有顶点及重心均被多边形包含。
func isRingInsidePolygon(ring orb.Ring, poly orb.Polygon) bool {
	if len(ring) < 3 {
		return false
	}
	for _, pt := range ring {
		if !planar.PolygonContains(poly, pt) {
			return false
		}
	}
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

// splitByClosedRing 用闭合环切割多边形，生成外环+孔洞结构和独立的内部多边形。
// 外部多边形保留原外环并添加孔洞，内部多边形以闭合环为外环。
func splitByClosedRing(poly orb.Polygon, ring orb.Ring, minArea float64) []orb.Polygon {
	outside := make(orb.Polygon, 0, len(poly)+1)
	outside = append(outside, poly[0])
	for i := 1; i < len(poly); i++ {
		outside = append(outside, poly[i])
	}

	holeRing := make(orb.Ring, len(ring))
	copy(holeRing, ring)
	if planar.Area(holeRing) > 0 {
		holeRing.Reverse()
	}
	outside = append(outside, holeRing)

	innerRing := make(orb.Ring, len(ring))
	copy(innerRing, ring)
	if planar.Area(innerRing) < 0 {
		innerRing.Reverse()
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

// normalizePolygonOrientation 规范化多边形环的方向。
// 外环强制为逆时针（CCW），孔洞强制为顺时针（CW），符合 GeoJSON 右手定则。
func normalizePolygonOrientation(poly orb.Polygon) orb.Polygon {
	if len(poly) == 0 {
		return poly
	}
	result := make(orb.Polygon, len(poly))
	result[0] = make(orb.Ring, len(poly[0]))
	copy(result[0], poly[0])
	if planar.Area(result[0]) < 0 {
		result[0].Reverse()
	}
	result[0] = utils.CloseRing(result[0])
	for i := 1; i < len(poly); i++ {
		result[i] = make(orb.Ring, len(poly[i]))
		copy(result[i], poly[i])
		if planar.Area(result[i]) > 0 {
			result[i].Reverse()
		}
		result[i] = utils.CloseRing(result[i])
	}
	return result
}

// minDistToRing 计算点到多边形环的最短距离。
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

// pointToSegmentDist 计算点到线段的最短欧几里得距离。
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

// lineLength 计算线段的累积长度（所有子段欧几里得距离之和）。
func lineLength(line orb.LineString) float64 {
	total := 0.0
	for i := 0; i < len(line)-1; i++ {
		total += math.Hypot(line[i+1][0]-line[i][0], line[i+1][1]-line[i][1])
	}
	return total
}

// splitByOpenLine 用开放线将多边形一分为二。
//
// 流程：
//  1. 找到线起止点在边界外环上的映射点
//  2. 将外环在映射点处拆为两段弧（arcAB 和 arcBA）
//  3. 提取线的正反两个方向段
//  4. 组合为两个新多边形
//  5. 将孔洞分配给包含它们的新多边形
func splitByOpenLine(poly orb.Polygon, line orb.LineString, minArea float64) []orb.Polygon {
	exterior := poly[0]

	if minDistToRing(line[0], exterior) > 0.01 || minDistToRing(line[len(line)-1], exterior) > 0.01 {
		return []orb.Polygon{poly}
	}

	firstPt, _ := nearestPointOnRing(line[0], poly[0])
	lastPt, _ := nearestPointOnRing(line[len(line)-1], poly[0])
	first := intersectionPt{pt: firstPt, lineIdx: 0, dist: 0}
	last := intersectionPt{pt: lastPt, lineIdx: len(line) - 2, dist: lineLength(line)}

	if utils.EqPoint(first.pt, last.pt, 0.001) {
		return []orb.Polygon{poly}
	}

	cutForward := extractLinePortion(line, first, last)
	cutBackward := reverseLineString(cutForward)

	arcAB, arcBA := splitRingAtPoints(poly[0], first.pt, last.pt)
	if arcAB == nil || arcBA == nil {
		return []orb.Polygon{poly}
	}

	poly1 := buildPolygonFromSegments(cutForward, arcBA)
	poly2 := buildPolygonFromSegments(cutBackward, arcAB)

	// 将孔洞分配给包含它们的新多边形
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

// nearestPointOnRing 查找点在多边形环上的最近点及对应距离。
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

// splitRingAtPoints 在两点处将环拆分为两段弧。
// 先将两点插入环的顶点序列，再从 idxA→idxB 和 idxB→idxA 分别提取两段弧。
func splitRingAtPoints(ring orb.Ring, pA, pB orb.Point) (orb.LineString, orb.LineString) {
	if len(ring) < 3 {
		return nil, nil
	}
	if minDistToRing(pA, ring) > 0.01 || minDistToRing(pB, ring) > 0.01 {
		return nil, nil
	}
	if utils.EqPoint(pA, pB, 0.001) {
		return nil, nil
	}

	n := len(ring)
	if utils.EqPoint(ring[0], ring[n-1], 1e-9) {
		n--
	}
	verts := ring[:n]

	expanded := insertPointsIntoVerts(verts, pA, pB)

	idxA := indexOfPoint(expanded, pA)
	idxB := indexOfPoint(expanded, pB)
	if idxA < 0 || idxB < 0 {
		return nil, nil
	}

	arcAB := extractArcForward(expanded, idxA, idxB)
	arcBA := extractArcForward(expanded, idxB, idxA)
	return arcAB, arcBA
}

// insertPointsIntoVerts 将两个点插入环的顶点序列中的正确位置。
// 先定位每个点所在的段（已有的顶点或投影到某段上），
// 若同一段上有多个插入点，按投影参数 t 排序。
func insertPointsIntoVerts(verts []orb.Point, pA, pB orb.Point) []orb.Point {
	n := len(verts)

	type locateResult struct {
		pt     orb.Point
		segIdx int
		t      float64
		onVtx  bool
	}

	locate := func(pt orb.Point) locateResult {
		for i, v := range verts {
			if utils.EqPoint(pt, v, 0.001) {
				return locateResult{pt, i, 0, true}
			}
		}
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

	segInserts := make([][]locateResult, n)
	for _, loc := range []locateResult{locA, locB} {
		if loc.onVtx {
			continue
		}
		segInserts[loc.segIdx] = append(segInserts[loc.segIdx], loc)
	}

	for i := range segInserts {
		if len(segInserts[i]) > 1 {
			sort.Slice(segInserts[i], func(a, b int) bool {
				return segInserts[i][a].t < segInserts[i][b].t
			})
		}
	}

	result := make([]orb.Point, 0, n+2)
	for i := 0; i < n; i++ {
		result = append(result, verts[i])
		for _, ins := range segInserts[i] {
			result = append(result, ins.pt)
		}
	}

	return result
}

// indexOfPoint 在顶点列表中查找指定点的索引（带容差）。
func indexOfPoint(verts []orb.Point, pt orb.Point) int {
	for i, v := range verts {
		if utils.EqPoint(v, pt, 0.001) {
			return i
		}
	}
	return -1
}

// extractArcForward 沿顶点列表正向提取从 from 到 to 的弧段。
// 支持跨列表尾部的环形提取。
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
	return pts
}

// extractLinePortion 提取等值线在两个交点之间的部分。
func extractLinePortion(line orb.LineString, first, last intersectionPt) orb.LineString {
	var pts orb.LineString
	pts = append(pts, first.pt)

	for li := first.lineIdx + 1; li <= last.lineIdx; li++ {
		if !utils.EqPoint(pts[len(pts)-1], line[li], 1e-9) {
			pts = append(pts, line[li])
		}
	}

	if !utils.EqPoint(pts[len(pts)-1], last.pt, 1e-9) {
		pts = append(pts, last.pt)
	}

	return pts
}

// reverseLineString 反转线段的方向。
func reverseLineString(ls orb.LineString) orb.LineString {
	n := len(ls)
	result := make(orb.LineString, n)
	for i := 0; i < n; i++ {
		result[i] = ls[n-1-i]
	}
	return result
}

// buildPolygonFromSegments 用两段弧拼接为一个多边形。
// 自动处理首尾去重和闭合，确保外环为逆时针方向。
func buildPolygonFromSegments(seg1, seg2 orb.LineString) orb.Polygon {
	var ring orb.Ring
	ring = append(ring, seg1...)
	if len(seg2) > 0 && len(ring) > 0 && utils.EqPoint(ring[len(ring)-1], seg2[0], 1e-9) {
		ring = append(ring, seg2[1:]...)
	} else {
		ring = append(ring, seg2...)
	}
	if len(ring) > 0 && !utils.EqPoint(ring[0], ring[len(ring)-1], 1e-9) {
		ring = append(ring, ring[0])
	}
	ring = cleanRing(ring)
	if planar.Area(ring) < 0 {
		ring.Reverse()
	}
	ring = utils.CloseRing(ring)
	if len(ring) < 4 {
		return orb.Polygon{ring}
	}
	return orb.Polygon{ring}
}

// cleanRing 移除环中连续的重复顶点（在容差范围内）。
func cleanRing(ring orb.Ring) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	result := orb.Ring{ring[0]}
	for i := 1; i < len(ring); i++ {
		if !utils.EqPoint(result[len(result)-1], ring[i], 1e-9) {
			result = append(result, ring[i])
		}
	}
	return result
}

// polyWithLevel 带 level 属性的多边形，用于并行计算 level 时的中间结构。
type polyWithLevel struct {
	poly  orb.Polygon
	level float64
}

// deduplicatePolygons 通过规范化后的字符串 key 去重多边形列表。
func deduplicatePolygons(polys []orb.Polygon) []orb.Polygon {
	if len(polys) <= 1 {
		return polys
	}
	seen := make(map[string]struct{})
	result := polys[:0]
	for _, poly := range polys {
		key := polygonKey(poly)
		if _, ok := seen[key]; ok {
			continue
		}
		seen[key] = struct{}{}
		result = append(result, poly)
	}
	return result
}

// polygonKey 为多边形生成唯一的字符串标识。
// 通过对每个环进行规范化（从最小顶点开始），确保旋转等价的多边形获得相同的 key。
func polygonKey(poly orb.Polygon) string {
	var b strings.Builder
	buf := make([]byte, 0, 24)
	for ri, ring := range poly {
		if ri > 0 {
			b.WriteByte('|')
		}
		normRing := normalizeRing(ring)
		for i, pt := range normRing {
			if i > 0 {
				b.WriteByte(';')
			}
			buf = strconv.AppendFloat(buf[:0], pt[0], 'f', 6, 64)
			b.Write(buf)
			b.WriteByte(',')
			buf = strconv.AppendFloat(buf[:0], pt[1], 'f', 6, 64)
			b.Write(buf)
		}
	}
	return b.String()
}

// normalizeRing 将环从最小的顶点开始旋转，确保等价环的规范化形式一致。
func normalizeRing(ring orb.Ring) orb.Ring {
	if len(ring) < 3 {
		return ring
	}
	n := len(ring) - 1
	if n < 2 {
		return ring
	}

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

// polygonSideOfLine 判断多边形相对于等值线的位置侧（left/right/inside/outside）。
// 首先检查是否与边界有共享边，若有则通过共享边方向判断；
// 若无共享边则通过点与线的叉积关系判断。
func polygonSideOfLine(poly orb.Polygon, entry contourEntry) (side string, hasSharedEdge bool) {
	if len(entry.coords) < 2 {
		return "", false
	}
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

	start, end, reversed, found := findSharedEdgeSegment(poly[0], entry.geom)

	if found {
		dx := end[0] - start[0]
		dy := end[1] - start[1]
		l := math.Hypot(dx, dy)
		if l < 1e-12 {
			return "right", true
		}
		mid := orb.Point{(start[0] + end[0]) / 2, (start[1] + end[1]) / 2}

		// 多点采样判断：沿左法向取多个偏移距离，多数投票。
		// 与 contour.go 的 sampleHighSide 多点采样思路一致，避免单点误判。
		offsets := []float64{0.0005, 0.001, 0.002, 0.004}
		leftCount := 0
		validCount := 0
		for _, eps := range offsets {
			leftPt := orb.Point{mid[0] - dy/l*eps, mid[1] + dx/l*eps}
			onLeft := planar.PolygonContains(poly, leftPt)
			if reversed {
				onLeft = !onLeft
			}
			validCount++
			if onLeft {
				leftCount++
			}
		}
		if leftCount > validCount/2 {
			return "left", true
		}
		return "right", true
	}

	return pointSideOfLine(repPt, entry.coords), true
}

// computePolygonLevel 为给定的多边形计算其等值线 level。
//
// 遍历所有 contourEntry，找到第一个与多边形共享边的等值线：
//   - 若多边形在该等值线的 highSide → 返回该 entry 的 level
//   - 否则返回上一个 level（prevLevel）
func computePolygonLevel(poly orb.Polygon, entries []contourEntry, opt *ContourOption) (level float64) {
	for _, entry := range entries {
		side, hasSharedEdge := polygonSideOfLine(poly, entry)
		if !hasSharedEdge {
			continue
		}
		if entry.highSide == side {
			return entry.level
		}
		return prevLevel(entry.level, opt)
	}
	return
}

// prevLevel 返回指定 level 在 LevelList 中的前一个级别。
// 若 LevelList 中无匹配，则用 ContourInterval 估算。
func prevLevel(level float64, opt *ContourOption) float64 {
	if len(opt.LevelList) > 0 {
		for i := 1; i < len(opt.LevelList); i++ {
			if opt.LevelList[i] == level {
				return opt.LevelList[i-1]
			}
		}
	}
	return math.Round((level-opt.ContourInterval)*1e4) / 1e4
}

// sharesEdgeWithExterior 判断多边形外环与等值线几何对象是否有共享边（总长度 > 0.01）。
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
			if (utils.EqPoint(ring[i], contourLS[j], 0.001) && utils.EqPoint(ring[i+1], contourLS[j+1], 0.001)) ||
				(utils.EqPoint(ring[i], contourLS[j+1], 0.001) && utils.EqPoint(ring[i+1], contourLS[j], 0.001)) {
				segLen := math.Hypot(ring[i+1][0]-ring[i][0], ring[i+1][1]-ring[i][1])
				totalMatch += segLen
			}
		}
	}
	return totalMatch > 0.01
}

// findSharedEdgeSegment 在多边形外环中查找与等值线几何对象共享的边段。
// 返回共享段的起点、终点、是否方向相反和是否找到。
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
			if utils.EqPoint(ring[i], contourLS[j], 0.001) && utils.EqPoint(ring[i+1], contourLS[j+1], 0.001) {
				return ring[i], ring[i+1], false, true
			} else if utils.EqPoint(ring[i], contourLS[j+1], 0.001) && utils.EqPoint(ring[i+1], contourLS[j], 0.001) {
				return ring[i], ring[i+1], true, true
			}
		}
	}
	return
}

// representativePoint 查找多边形内部的代表点。
// 优先使用多边形重心，若重心在多边形外部，则用第一条边的中点向内侧偏移。
func representativePoint(poly orb.Polygon) orb.Point {
	if len(poly) == 0 || len(poly[0]) < 4 {
		return orb.Point{0, 0}
	}

	cx, cy := 0.0, 0.0
	n := len(poly[0]) - 1
	for i := 0; i < n; i++ {
		cx += poly[0][i][0]
		cy += poly[0][i][1]
	}
	cx /= float64(n)
	cy /= float64(n)

	center := orb.Point{cx, cy}
	if planar.PolygonContains(poly, center) {
		return center
	}

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

// pointSideOfLine 使用叉积法判断点在线的哪一侧。
// cross > 0 → left, cross < 0 → right。
// pointSideOfLine 多点采样判断点 pt 在多段线的哪一侧（多数投票）。
// 沿多段线取 numSamples 条均匀分布的线段，对每条线段用叉积判断 pt 在其哪一侧，
// 过半左侧则返回 "left"，否则 "right"。
func pointSideOfLine(pt orb.Point, coords []*orb.Point) string {
	n := len(coords)
	if n < 2 {
		return "right"
	}

	numSamples := 5
	if n-1 < numSamples {
		numSamples = n - 1
	}
	if numSamples < 1 {
		numSamples = 1
	}

	stepSize := float64(n-1) / float64(numSamples)
	leftCount := 0
	validCount := 0

	for i := 0; i < numSamples; i++ {
		idx := int(float64(i) * stepSize)
		if idx >= n-1 {
			idx = n - 2
		}
		sx, sy := coords[idx][0], coords[idx][1]
		ex, ey := coords[idx+1][0], coords[idx+1][1]
		dx := ex - sx
		dy := ey - sy
		if math.Abs(dx) < 1e-9 && math.Abs(dy) < 1e-9 {
			continue
		}
		px := pt[0] - sx
		py := pt[1] - sy
		cross := dx*py - dy*px
		validCount++
		if cross > 0 {
			leftCount++
		}
	}

	if validCount == 0 || leftCount > validCount/2 {
		return "left"
	}
	return "right"
}

// polygonToFeature 将多边形转换为 GeoJSON Feature。
// 包括规范化环方向、四舍五入坐标（保留 6 位小数）和添加属性。
func polygonToFeature(poly orb.Polygon, level float64, opt *ContourOption) *geojson.Feature {
	if len(poly) == 0 {
		return nil
	}

	poly = normalizePolygonOrientation(poly)

	exterior := poly[0]
	if len(exterior) < 4 {
		return nil
	}

	resultPoly := make(orb.Polygon, 0, len(poly))
	resultPoly = append(resultPoly, roundRing(exterior))
	for i := 1; i < len(poly); i++ {
		resultPoly = append(resultPoly, roundRing(poly[i]))
	}

	feature := geojson.NewFeature(resultPoly)
	feature.Properties = maputil.Merge(map[string]interface{}{
		"type":  opt.ContourType,
		"level": math.Round(level*1e4) / 1e4,
	}, opt.Extra)

	return feature
}

// roundRing 将环的所有顶点四舍五入到 6 位小数。
func roundRing(ring orb.Ring) orb.Ring {
	result := make(orb.Ring, len(ring))
	for i, pt := range ring {
		result[i] = orb.Point{math.Round(pt[0]*1e6) / 1e6, math.Round(pt[1]*1e6) / 1e6}
	}
	return result
}
