// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了等值线生成的算法实现，包括网格插值、等值线提取、几何处理等。
package kriging_contour

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/duke-git/lancet/v2/maputil"
	"github.com/inkbamboo/kriging-contour/internal/utils"
	"github.com/panjf2000/ants/v2"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	"github.com/paulmach/orb/planar"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

// plotMutex 保护 gonum/plot 的渲染操作，因为 plotter.NewContour 不是并发安全的。
var plotMutex sync.Mutex

// tolerance 线段端点合并的容差距离。用于 mergeConnectedLines 中判断两条线段是否端点重合。
var tolerance = 0.01

// extractCoords 将 vg.Path 中的点坐标从画布空间转换回数据空间。
// canvasSize 为画布尺寸，xMin/xMax/yMin/yMax 为数据范围。
// 返回经四舍五入（保留 6 位小数）的地理坐标点序列。
func extractCoords(p vg.Path, xMin, xMax, yMin, yMax, canvasSize float64) []*orb.Point {
	var pts []*orb.Point
	for _, comp := range p {
		if comp.Type != vg.MoveComp && comp.Type != vg.LineComp {
			continue
		}
		x := float64(comp.Pos.X)/canvasSize*(xMax-xMin) + xMin
		y := float64(comp.Pos.Y)/canvasSize*(yMax-yMin) + yMin
		pts = append(pts, &orb.Point{math.Round(x*1e6) / 1e6, math.Round(y*1e6) / 1e6})
	}
	return pts
}

// GenerateLines 从散点数据生成等值线。这是等值线生成的主入口函数。
//
// 完整流程：
//  1. 输出数据统计信息
//  2. 校验参数和数据
//  3. 执行克里金插值生成网格（GenerateGrid）
//  4. 从网格中提取等值线（extractLinesFromGrid）
//
// 参数：
//   - points: 离散采样点数据
//   - boundary: 边界多边形，用于裁剪等值线和限制插值范围
//   - opt: 等值线生成配置参数
//
// 返回 GeoJSON Feature 切片，每个 Feature 代表一条等值线。
func GenerateLines(points []*Point, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature {
	stats := ComputeDataStats(points)
	fmt.Printf("[INFO] 原始数据点: %d, 有效点: %d, X范围: [%.4f, %.4f], Y范围: [%.4f, %.4f], Z范围: [%.4f, %.4f], NaN: %d, Inf: %d, 重复: %d\n",
		stats.OriginalCount, stats.CleanedCount,
		stats.MinX, stats.MaxX, stats.MinY, stats.MaxY,
		stats.MinZ, stats.MaxZ, stats.NanCount, stats.InfCount, stats.DupCount)
	var err error

	if points, err = CheckParams(points, boundary, opt); err != nil {
		return nil
	}
	var grid *Grid
	if grid, err = GenerateGrid(points, boundary, opt.Resolution, opt.LevelList); err != nil || grid == nil {
		fmt.Printf("[ERROR] 网格生成失败: %v\n", err)
		return nil
	}

	return extractLinesFromGrid(grid, boundary, opt)
}

// CheckParams 串联调用 ValidatePoints、ValidateContourOption 和 ValidateBoundary，
// 对数据、参数和边界进行完整校验。返回清洗后的数据点。
func CheckParams(points []*Point, boundary orb.Polygon, opt *ContourOption) (cleanPoints []*Point, err error) {
	if len(points) == 0 {
		fmt.Printf("[ERROR] 数据为空\n")
		return nil, fmt.Errorf("数据为空")
	}
	if opt == nil {
		fmt.Printf("[ERROR] 参数为空\n")
		return nil, fmt.Errorf("参数为空")
	}
	if cleanPoints, err = ValidatePoints(points, 3); err != nil {
		fmt.Printf("[ERROR] 数据校验失败: %v\n", err)
		return
	}
	if err = ValidateContourOption(opt); err != nil {
		fmt.Printf("[ERROR] 参数校验失败: %v\n", err)
		return
	}
	if err = ValidateBoundary(boundary); err != nil {
		fmt.Printf("[ERROR] 边界校验失败: %v\n", err)
		return
	}
	// 根据清洗后的数据补全 LevelList（若未显式指定）
	if len(opt.LevelList) == 0 {
		opt.LevelList = GetLevelList(cleanPoints, opt)
	}
	return
}

// GetLevelList 根据散点数据和配置参数计算等值线级别列表。
//
// 从 points 中提取有效的 Z 值范围（过滤 nil/NaN/Inf），
// 然后调用 opt.GetLevelList(minZ, maxZ) 生成级别列表。
func GetLevelList(points []*Point, opt *ContourOption) []float64 {
	if len(opt.LevelList) > 0 {
		return opt.LevelList
	}

	minZ, maxZ := math.Inf(1), math.Inf(-1)
	firstValid := true

	for _, p := range points {
		if p == nil || math.IsNaN(p.Z) || math.IsInf(p.Z, 0) {
			continue
		}
		if firstValid {
			minZ, maxZ = p.Z, p.Z
			firstValid = false
		} else {
			if p.Z < minZ {
				minZ = p.Z
			}
			if p.Z > maxZ {
				maxZ = p.Z
			}
		}
	}
	start := opt.ContourStart
	end := opt.ContourEnd

	if start == 0 && end == 0 {
		start = math.Floor(minZ/opt.ContourInterval) * opt.ContourInterval
		if start >= minZ {
			start -= opt.ContourInterval
		}
		end = math.Ceil(maxZ/opt.ContourInterval) * opt.ContourInterval
		if end <= maxZ {
			end += opt.ContourInterval
		}
	}

	if start > end {
		start, end = end, start
	}

	var list []float64
	for lv := start; lv <= end+1e-10; lv += opt.ContourInterval {
		list = append(list, math.Round(lv*1e4)/1e4)
	}
	return list
}

// GenerateLinesByGrid 基于已有的插值网格直接提取等值线。
// 适用于需要自定义插值参数或复用已有网格的场景。
func GenerateLinesByGrid(grid *Grid, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature {
	if grid == nil {
		fmt.Printf("[ERROR] 网格为空\n")
		return nil
	}

	return extractLinesFromGrid(grid, boundary, opt)
}

// extractLinesFromGrid 从插值网格中提取所有等值线。
//
// 处理步骤：
//  1. 确保边界多边形闭合
//  2. 使用 goroutine 池并行处理每个 level
//  3. 对每条等值线进行合并、裁剪、定向
//  4. 附加边界 Feature
//  5. 可选渲染 PNG 图片
func extractLinesFromGrid(grid *Grid, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature {
	fmt.Printf("[INFO] 等值线级别数: %d\n", len(opt.LevelList))
	closeBoundaryPolygon(boundary)
	extractStart := time.Now()
	fmt.Printf("[INFO] 步骤开始: 等值线提取\n")
	if len(opt.LevelList) == 0 {
		fmt.Printf("[WARN] 没有有效的等值线级别\n")
		fmt.Printf("[WARN] 步骤完成(有警告): 等值线提取 (耗时: %v)\n", time.Since(extractStart).Round(time.Millisecond))
		return nil
	}

	rows, cols := grid.Data.Dims()
	if rows < 2 || cols < 2 {
		fmt.Printf("[ERROR] 网格尺寸过小: %d × %d\n", rows, cols)
		fmt.Printf("[ERROR] 步骤失败: 等值线提取 (耗时: %v) - 网格尺寸过小\n", time.Since(extractStart).Round(time.Millisecond))
		return nil
	}
	const canvasSize = 2000.0

	var features []*geojson.Feature

	var wg sync.WaitGroup
	var lock sync.Mutex
	baseProperties := maputil.Merge(map[string]interface{}{
		"type": opt.ContourType,
	}, opt.Extra)
	// 使用 goroutine 池并行处理每个 level 的等值线提取
	p, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(body interface{}) {
		defer wg.Done()
		levelFeatures := generateLevelLines(grid, body.(float64), canvasSize, boundary, baseProperties)
		lock.Lock()
		features = append(features, levelFeatures...)
		defer lock.Unlock()
	})
	defer p.Release()
	for _, level := range opt.LevelList {
		wg.Add(1)
		_ = p.Invoke(level)
	}
	wg.Wait()

	fmt.Printf("[INFO] 步骤完成: 等值线提取 (耗时: %v, features: %d)\n", time.Since(extractStart).Round(time.Millisecond), len(features))
	fmt.Printf("[INFO] 等值线总数: %d\n", len(features))

	// 附加边界 Feature
	boundaryFeature := geojson.NewFeature(boundary)
	boundaryFeature.Properties = maputil.Merge(map[string]interface{}{"type": "boundary"}, opt.Extra)
	features = append(features, boundaryFeature)

	if opt.ImagePath != "" {
		imgStart := time.Now()
		fmt.Printf("[INFO] 步骤开始: 渲染等值线图片%+v\n", opt.LevelList)
		renderContourImage(grid, opt.LevelList, opt.ImagePath)
		fmt.Printf("[INFO] 步骤完成: 渲染等值线图片 (耗时: %v)\n", time.Since(imgStart).Round(time.Millisecond))
	}

	return features
}

// closeBoundaryPolygon 确保多边形的所有环都是闭合的（首尾点相同）。
func closeBoundaryPolygon(poly orb.Polygon) {
	for ri, ring := range poly {
		if len(ring) > 0 && ring[0] != ring[len(ring)-1] {
			poly[ri] = append(ring, ring[0])
		}
	}
}

// generateLevelLines 为单个等值线 level 生成 GeoJSON Feature 列表。
//
// 处理流程：
//  1. 使用 gonum/plot 的 NewContour（Marching Squares）在自定义 canvas 上绘制
//  2. 从 vg.Path 提取坐标并转换回地理空间
//  3. 合并端点相接的线段
//  4. 用边界多边形裁剪
//  5. 过滤过短的线段
//  6. 统一线段方向（resolveLineOrientationCG）
//  7. 交换 X/Y 坐标以匹配 orb 库的 (lon, lat) 约定
func generateLevelLines(grid *Grid, level float64, canvasSize font.Length, boundary orb.Polygon, baseProperties map[string]interface{}) []*geojson.Feature {
	var features []*geojson.Feature
	c := plotter.NewContour(grid, []float64{level}, nil)
	xMin, xMax, yMin, yMax := c.DataRange()

	rc := &canvas{}
	dc := draw.Canvas{
		Canvas: rc,
		Rectangle: vg.Rectangle{
			Min: vg.Point{},
			Max: vg.Point{X: canvasSize, Y: canvasSize},
		},
	}
	p := plot.New()
	p.X.Min, p.X.Max = xMin, xMax
	p.Y.Min, p.Y.Max = yMin, yMax
	p.Add(c)
	// 串行化 gonum/plot 的渲染，因其内部不是并发安全的
	plotMutex.Lock()
	c.Plot(dc, p)
	plotMutex.Unlock()

	var allCoords [][]*orb.Point
	for _, path := range rc.getPaths() {
		coords := extractCoords(path, xMin, xMax, yMin, yMax, float64(canvasSize))
		if len(coords) < 2 {
			continue
		}
		allCoords = append(allCoords, coords)
	}

	allCoords = clipLinesToPolygon(allCoords, boundary)

	merged := mergeConnectedLines(allCoords)
	// 过滤长度小于边界周长 1/100000 的线段
	minLength := planar.Length(orb.LineString(boundary[0])) * 0.00001
	for _, coords := range merged {
		if len(coords) < 2 || planar.Length(orb.LineString(toOrbPoints(coords))) < minLength {
			continue
		}

		isClosed, highSide := resolveLineOrientationCG(grid, coords, level, boundary)

		// 坐标交换 (lon,lat)→(lat,lon)，匹配 orb 输出的 (lat,lon) 约定
		// 交换 X/Y 后几何方向反转：(lon,lat) 的 CCW 在 (lat,lon) 变为 CW
		for i := range coords {
			coords[i][0], coords[i][1] = coords[i][1], coords[i][0]
		}

		// 补偿坐标交换的影响：
		//   封闭线：CCW→CW，反转恢复 CCW
		//   开放线：left↔right 翻转
		if isClosed {
			reversePoints(coords)
		} else {
			if highSide == "left" {
				highSide = "right"
			} else {
				highSide = "left"
			}
		}

		ls := orb.LineString(toOrbPoints(coords))
		feature := geojson.NewFeature(ls)
		feature.Properties = maputil.Merge(map[string]interface{}{}, baseProperties)
		feature.Properties["level"] = math.Round(level*1e4) / 1e4
		feature.Properties["is_closed"] = isClosed
		feature.Properties["high_value_side"] = highSide
		features = append(features, feature)
	}
	return features
}

// toOrbPoints 将 []*orb.Point 转换为 []orb.Point。
func toOrbPoints(pts []*orb.Point) []orb.Point {
	result := make([]orb.Point, len(pts))
	for i, p := range pts {
		result[i] = *p
	}
	return result
}

// fromOrbPoints 将 []orb.Point 转换为 []*orb.Point。
func fromOrbPoints(pts []orb.Point) []*orb.Point {
	result := make([]*orb.Point, len(pts))
	for i := range pts {
		result[i] = &orb.Point{pts[i][0], pts[i][1]}
	}
	return result
}

// resolveLineOrientationCG 判断等值线的几何闭合状态及高值侧方向。
//
// 统一采样逻辑（封闭/非封闭共用）：
//
//	取多段线中间有向线段的中点，沿左法向（行进方向逆时针90°）偏移微小距离，
//	采样该点的网格插值。比较采样值与 level 的大小即可得出高值侧。
//
// 封闭线：
//
//	先确保多边形为 CCW（逆时针 / 左旋），此时有向线段的左侧 = 多边形内部。
//	左法向采样值 > level → 高值在内部 → "inside"
//	左法向采样值 ≤ level → 高值在外部 → "outside"
//
// 非封闭线：
//
//	开线将边界切割成两个子多边形。方向需保证从起点沿开线到终点，
//	再沿边界左旋（CCW）回到起点所围成的封闭多边形为 CCW。
//	若该多边形面积 < 0（CW），则反转开线方向。
//	左法向采样值 > level → 高值在左侧 → "left"
//	左法向采样值 ≤ level → 高值在右侧 → "right"
//
// 注意：调用方在 (lon,lat)→(lat,lon) 坐标交换后需对 highSide 做翻转补偿。
func resolveLineOrientationCG(grid *Grid, coords []*orb.Point, level float64, boundary orb.Polygon) (isClosed bool, highSide string) {
	closed := isLineClosed(coords, 0.01)

	if closed {
		// 封闭线：确保 CCW（左旋），使有向线段的左侧指向多边形内部
		if planar.Area(orb.Ring(toOrbPoints(coords))) < 0 {
			reversePoints(coords)
		}
	} else if len(boundary) > 0 && len(boundary[0]) >= 3 {
		// 非封闭线：通过构造 (开线 + 边界段) 的封闭多边形来判断方向
		// 开线方向需使得 起点→终点→边界左旋→起点 围成的多边形为 CCW
		ring := boundary[0]
		startPt := *coords[0]
		endPt := *coords[len(coords)-1]

		// 找起终点在边界环上的最近段
		si := closestRingSegment(startPt, ring)
		ei := closestRingSegment(endPt, ring)

		// 沿边界环从终点左旋（CCW）回到起点，构建边界段
		var ringPts []orb.Point
		if si == ei {
			// 起终点在同一段上，直接用开线本身构成闭合环来判方向
			// （开线两端 + 沿边界同一段回到起点）
			ringPts = append(ringPts, endPt)
			for j := ei; j != si; j = (j + 1) % (len(ring) - 1) {
				ringPts = append(ringPts, ring[j])
			}
			ringPts = append(ringPts, startPt)
		} else {
			// 从终点段末 → 沿环到起点段首
			for j := ei; j != si; j = (j + 1) % (len(ring) - 1) {
				ringPts = append(ringPts, ring[j])
			}
			ringPts = append(ringPts, ring[si])
		}

		// 构建封闭多边形：开线坐标 + 边界段
		polyPts := make([]orb.Point, 0, len(coords)+len(ringPts))
		for _, p := range coords {
			polyPts = append(polyPts, *p)
		}
		polyPts = append(polyPts, ringPts...)
		// 确保闭合
		if !utils.EqPoint(polyPts[0], polyPts[len(polyPts)-1], 1e-9) {
			polyPts = append(polyPts, polyPts[0])
		}

		if planar.Area(orb.Ring(polyPts)) < 0 {
			reversePoints(coords)
		}
	}

	// 多点采样：在多段线的多个位置沿左法向采样，多数投票决定高值侧。
	// 避免单点采样因局部异常导致误判。
	highSide = sampleHighSide(grid, coords, level, closed)
	return closed, highSide
}

// sampleHighSide 采样网格值，多数投票判断高值侧。
//
// 封闭线：优先用质心采样（最稳定），若质心在多边形外（极度凹多边形），
//
//	则退化为沿多段线左法向多点采样（CCW 封闭线左侧=内部）。
//
// 非封闭线：沿多段线多点左法向采样，过半则为 "left"，否则 "right"。
func sampleHighSide(grid *Grid, coords []*orb.Point, level float64, closed bool) string {
	n := len(coords)
	if n < 2 {
		if closed {
			return "inside"
		}
		return "left"
	}

	if closed {
		// 质心采样（最可靠的方式）
		cx, cy := 0.0, 0.0
		for i := 0; i < n; i++ {
			cx += coords[i][0]
			cy += coords[i][1]
		}
		cx /= float64(n)
		cy /= float64(n)

		ring := orb.Ring(toOrbPoints(coords))
		if planar.PolygonContains(orb.Polygon{ring}, orb.Point{cx, cy}) {
			if grid.Sample(cx, cy) > level {
				return "inside"
			}
			return "outside"
		}
		// 质心在外（极度凹多边形），退化为边缘采样
	}

	// 多点边缘采样：取 numSamples 条均匀线段，中点沿左法向偏移采样，多数投票
	numSamples := 5
	if n-1 < numSamples {
		numSamples = n - 1
	}
	if numSamples < 1 {
		numSamples = 1
	}

	stepSize := float64(n-1) / float64(numSamples)
	countHigh := 0
	countValid := 0

	for i := 0; i < numSamples; i++ {
		idx := int(float64(i) * stepSize)
		if idx >= n-1 {
			idx = n - 2
		}
		a, b := coords[idx], coords[idx+1]
		dx := b[0] - a[0] // Δlon
		dy := b[1] - a[1] // Δlat

		segLen := math.Hypot(dx, dy)
		if segLen < 1e-9 {
			continue
		}

		// 中点沿左法向 (-dy, dx) 偏移后采样
		mx := (a[0] + b[0]) / 2
		my := (a[1] + b[1]) / 2
		step := 0.1
		val := grid.Sample(mx-dy*step, my+dx*step)

		countValid++
		if val > level {
			countHigh++
		}
	}

	if countValid == 0 {
		if closed {
			return "inside"
		}
		return "left"
	}

	if countHigh > countValid/2 {
		if closed {
			return "inside"
		}
		return "left"
	}
	if closed {
		return "outside"
	}
	return "right"
}

// closestRingSegment 找到点 pt 在环 ring 上最近的段索引。
// 返回段起点在环中的索引。
func closestRingSegment(pt orb.Point, ring orb.Ring) int {
	minDist := math.MaxFloat64
	minIdx := 0
	for i := 0; i < len(ring)-1; i++ {
		cp := closestPointOnSegment(pt, ring[i], ring[i+1])
		d := (cp[0]-pt[0])*(cp[0]-pt[0]) + (cp[1]-pt[1])*(cp[1]-pt[1])
		if d < minDist {
			minDist = d
			minIdx = i
		}
	}
	return minIdx
}

// reversePoints 原地反转点序列。
func reversePoints(pts []*orb.Point) {
	for i, j := 0, len(pts)-1; i < j; i, j = i+1, j-1 {
		pts[i], pts[j] = pts[j], pts[i]
	}
}

// mergeConnectedLines 合并端点相接的线段。
//
// 使用空间哈希（cell grid）加速邻近搜索：
//  1. 将每条线的首尾端点映射到空间格网
//  2. 在同一格网及相邻格网内寻找距离在 tolerance 内的端点对
//  3. 使用并查集（Union-Find）将属于同一线段组的线索引合并
//  4. 对每个连通分量调用 rebuildChain 重建有序线段链
func mergeConnectedLines(lines [][]*orb.Point) [][]*orb.Point {
	if len(lines) <= 1 {
		return lines
	}

	n := len(lines)
	cellSize := tolerance

	type epEntry struct {
		idx    int
		isHead bool
	}

	epMap := make(map[string][]epEntry)
	addToCells := func(pt *orb.Point, entry epEntry) {
		cx := int(math.Floor(pt[0] / cellSize))
		cy := int(math.Floor(pt[1] / cellSize))
		for dx := -1; dx <= 1; dx++ {
			for dy := -1; dy <= 1; dy++ {
				k := fmt.Sprintf("%d,%d", cx+dx, cy+dy)
				epMap[k] = append(epMap[k], entry)
			}
		}
	}

	for i, line := range lines {
		if len(line) < 2 {
			continue
		}
		addToCells(line[0], epEntry{i, true})
		addToCells(line[len(line)-1], epEntry{i, false})
	}

	// 并查集
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	find := func(x int) int {
		root := x
		for parent[root] != root {
			root = parent[root]
		}
		for parent[x] != root {
			parent[x], x = root, parent[x]
		}
		return root
	}

	seen := make(map[[2]int]bool)
	for _, entries := range epMap {
		for i := 0; i < len(entries); i++ {
			for j := i + 1; j < len(entries); j++ {
				ea, eb := entries[i], entries[j]
				if ea.idx == eb.idx {
					continue
				}
				a, b := ea.idx, eb.idx
				if a > b {
					a, b = b, a
				}
				pair := [2]int{a, b}
				if seen[pair] {
					continue
				}
				seen[pair] = true

				var pa, pb *orb.Point
				if ea.isHead {
					pa = lines[ea.idx][0]
				} else {
					pa = lines[ea.idx][len(lines[ea.idx])-1]
				}
				if eb.isHead {
					pb = lines[eb.idx][0]
				} else {
					pb = lines[eb.idx][len(lines[eb.idx])-1]
				}
				if utils.EqPoint(*pa, *pb, tolerance) {
					ra, rb := find(ea.idx), find(eb.idx)
					if ra != rb {
						parent[ra] = rb
					}
				}
			}
		}
	}

	compMap := make(map[int][]int)
	for i := 0; i < n; i++ {
		if len(lines[i]) < 2 {
			continue
		}
		root := find(i)
		compMap[root] = append(compMap[root], i)
	}

	var result [][]*orb.Point
	for _, indices := range compMap {
		merged := rebuildChain(lines, indices)
		result = append(result, merged...)
	}

	return result
}

// rebuildChain 将属于同一组的线段片段重建为有序的线段链。
//
// 首先匹配各片段的头尾连接关系，然后从度数为 1 的端点（或闭合环的任意点）出发，
// 按连接关系依次拼接，最终输出一条或多条连续线段。
func rebuildChain(lines [][]*orb.Point, indices []int) [][]*orb.Point {
	if len(indices) == 1 {
		return [][]*orb.Point{lines[indices[0]]}
	}

	m := len(indices)
	idxToPos := make(map[int]int, m)
	posToIdx := make([]int, m)
	for pos, idx := range indices {
		idxToPos[idx] = pos
		posToIdx[pos] = idx
	}

	type conn struct {
		toPos  int
		toHead bool
	}
	headConn := make([]*conn, m)
	tailConn := make([]*conn, m)

	for a := 0; a < m; a++ {
		la := lines[posToIdx[a]]
		for b := a + 1; b < m; b++ {
			lb := lines[posToIdx[b]]

			if utils.EqPoint(*la[0], *lb[0], tolerance) {
				headConn[a] = &conn{b, true}
				headConn[b] = &conn{a, true}
			}
			if utils.EqPoint(*la[0], *lb[len(lb)-1], tolerance) {
				headConn[a] = &conn{b, false}
				tailConn[b] = &conn{a, true}
			}
			if utils.EqPoint(*la[len(la)-1], *lb[0], tolerance) {
				tailConn[a] = &conn{b, true}
				headConn[b] = &conn{a, false}
			}
			if utils.EqPoint(*la[len(la)-1], *lb[len(lb)-1], tolerance) {
				tailConn[a] = &conn{b, false}
				tailConn[b] = &conn{a, false}
			}
		}
	}

	degree := make([]int, m)
	for i := 0; i < m; i++ {
		if headConn[i] != nil {
			degree[i]++
		}
		if tailConn[i] != nil {
			degree[i]++
		}
	}

	visited := make([]bool, m)
	var result [][]*orb.Point

	for {
		// 优先从度数为 1 的端点开始（开放链的起点）
		start := -1
		startFromHead := true
		for i := 0; i < m; i++ {
			if visited[i] {
				continue
			}
			if degree[i] == 1 {
				start = i
				startFromHead = headConn[i] == nil
				break
			}
		}
		// 若无度数为 1 的端点，任选一个未访问的片段（闭合环）
		if start == -1 {
			for i := 0; i < m; i++ {
				if !visited[i] {
					start = i
					startFromHead = true
					break
				}
			}
			if start == -1 {
				break
			}
		}

		var chain []*orb.Point
		cur := start
		forward := startFromHead

		for {
			visited[cur] = true
			line := lines[posToIdx[cur]]

			if forward {
				chain = append(chain, line...)
			} else {
				for i := len(line) - 1; i >= 0; i-- {
					chain = append(chain, line[i])
				}
			}

			var next *conn
			if forward {
				next = tailConn[cur]
			} else {
				next = headConn[cur]
			}
			if next == nil || visited[next.toPos] {
				break
			}

			forward = next.toHead
			cur = next.toPos
		}

		if len(chain) >= 2 {
			result = append(result, chain)
		}
	}

	return result
}

// clipLinesToPolygon 用边界多边形裁剪所有等值线段。
// 先将边界多边形做坐标交换（(lat,lon)→(lon,lat)），再逐条裁剪。
func clipLinesToPolygon(lines [][]*orb.Point, boundary orb.Polygon) [][]*orb.Point {
	normPoly := make(orb.Polygon, len(boundary))
	for ri, ring := range boundary {
		normRing := make(orb.Ring, len(ring))
		for i, p := range ring {
			normRing[i] = orb.Point{p[1], p[0]}
		}
		normPoly[ri] = normRing
	}

	var result [][]*orb.Point
	for _, line := range lines {
		result = append(result, clipLineToPolygon(line, normPoly)...)
	}
	return result
}

// clipLineToPolygon 使用 Sutherland–Hodgman 风格算法裁剪单条线段。
//
// 对线段上的每对连续点，根据其在多边形内/外的状态分 4 种情况处理：
//   - 内→内：保留终点
//   - 内→外：找交点，添加到当前段，结束该段
//   - 外→内：找交点，开始新段
//   - 外→外：可能跨越整个多边形，找两个交点生成穿越段
func clipLineToPolygon(line []*orb.Point, boundary orb.Polygon) [][]*orb.Point {
	if len(line) < 2 {
		return nil
	}
	ring := boundary[0]
	var result [][]*orb.Point
	var current []*orb.Point

	prevIn := planar.PolygonContains(boundary, *line[0])
	if prevIn {
		current = append(current, line[0])
	}
	for i := 1; i < len(line); i++ {
		currIn := planar.PolygonContains(boundary, *line[i])
		if prevIn && currIn {
			current = append(current, line[i])
		} else if prevIn && !currIn {
			if pt, ok := lineRingFirstIntersection(*line[i-1], *line[i], ring); ok {
				current = append(current, &pt)
			}
			if len(current) >= 2 {
				result = append(result, current)
			}
			current = nil
		} else if !prevIn && currIn {
			if pt, ok := lineRingFirstIntersection(*line[i-1], *line[i], ring); ok {
				current = []*orb.Point{&pt, line[i]}
			} else {
				current = []*orb.Point{line[i]}
			}
		} else {
			crossPts := linePolygonCrossings(*line[i-1], *line[i], ring)
			if len(crossPts) == 2 {
				result = append(result, []*orb.Point{crossPts[0], crossPts[1]})
			}
		}

		prevIn = currIn
	}

	if len(current) >= 2 {
		result = append(result, current)
	}

	return result
}

// lineRingFirstIntersection 找线段与多边形环的最近交点（从线段起点 a 方向出发）。
// 返回距离 a 最近的交点和是否找到。
func lineRingFirstIntersection(a, b orb.Point, ring orb.Ring) (orb.Point, bool) {
	var bestPt orb.Point
	bestDist := math.MaxFloat64
	found := false

	for i := 0; i < len(ring)-1; i++ {
		pt, ok := segmentIntersection(a, b, ring[i], ring[i+1])
		if ok {
			d := (pt[0]-a[0])*(pt[0]-a[0]) + (pt[1]-a[1])*(pt[1]-a[1])
			if d < bestDist {
				bestDist = d
				bestPt = pt
				found = true
			}
		}
	}

	return bestPt, found
}

// closestPointOnSegment 计算点 p 到线段 ab 上的最近点（垂足或端点）。
func closestPointOnSegment(p, a, b orb.Point) orb.Point {
	dx, dy := b[0]-a[0], b[1]-a[1]
	lenSq := dx*dx + dy*dy
	if lenSq < 1e-15 {
		return a
	}
	t := ((p[0]-a[0])*dx + (p[1]-a[1])*dy) / lenSq
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	return orb.Point{a[0] + t*dx, a[1] + t*dy}
}

// linePolygonCrossings 找线段 ab 与多边形环的所有交点，返回最近的 2 个不同交点。
// 用于处理线段完全位于多边形外部但跨越整个多边形的情况。
func linePolygonCrossings(a, b orb.Point, ring orb.Ring) []*orb.Point {
	type ptDist struct {
		pt   *orb.Point
		dist float64
	}
	var crossings []ptDist

	for i := 0; i < len(ring)-1; i++ {
		pt, ok := segmentIntersection(a, b, ring[i], ring[i+1])
		if ok {
			d := (pt[0]-a[0])*(pt[0]-a[0]) + (pt[1]-a[1])*(pt[1]-a[1])
			crossings = append(crossings, ptDist{&orb.Point{pt[0], pt[1]}, d})
		}
	}

	if len(crossings) < 2 {
		return nil
	}

	sort.Slice(crossings, func(i, j int) bool {
		return crossings[i].dist < crossings[j].dist
	})

	var result []*orb.Point
	for i := 0; i < len(crossings); i++ {
		if len(result) == 0 {
			result = append(result, crossings[i].pt)
		} else {
			last := result[len(result)-1]
			d2 := (crossings[i].pt[0]-last[0])*(crossings[i].pt[0]-last[0]) +
				(crossings[i].pt[1]-last[1])*(crossings[i].pt[1]-last[1])
			if d2 > 1e-6 {
				result = append(result, crossings[i].pt)
			}
		}
		if len(result) >= 2 {
			break
		}
	}

	if len(result) < 2 {
		return nil
	}
	return result
}

// segmentIntersection 使用参数法计算两条线段的交点。
// 返回交点和是否相交（参数 t, u 均在 [−1e-12, 1+1e-12] 范围内）。
func segmentIntersection(p1, p2, p3, p4 orb.Point) (orb.Point, bool) {
	denom := (p1[0]-p2[0])*(p3[1]-p4[1]) - (p1[1]-p2[1])*(p3[0]-p4[0])
	if math.Abs(denom) < 1e-15 {
		return orb.Point{}, false
	}

	t := ((p1[0]-p3[0])*(p3[1]-p4[1]) - (p1[1]-p3[1])*(p3[0]-p4[0])) / denom
	u := -((p1[0]-p2[0])*(p1[1]-p3[1]) - (p1[1]-p2[1])*(p1[0]-p3[0])) / denom

	if t >= -1e-12 && t <= 1+1e-12 && u >= -1e-12 && u <= 1+1e-12 {
		return orb.Point{p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1])}, true
	}
	return orb.Point{}, false
}

// isLineClosed 判断线段是否闭合（首尾点距离在容差范围内且点数 >= 3）。
func isLineClosed(coords []*orb.Point, tol float64) bool {
	if len(coords) < 3 {
		return false
	}
	return utils.EqPoint(*coords[0], *coords[len(coords)-1], tol)
}
