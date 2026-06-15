// Package kriging_contour 基于克里金插值（Ordinary Kriging）从离散空间数据点生成等值线（LineString）和等值面（Polygon）。
//
// 主要功能：
//   - generateGrid: 使用克里金插值将离散点转换为规则网格
//   - GenerateLines: 从网格中提取等值线（marching squares 算法），含裁剪、合并、方向统一、短线过滤
//   - GenerateFaces: 基于等值线进行遍历分割生成等值面
//
// 等值线方向统一遵循右手定则：沿线条走向，高值区在左侧。
// 坐标系统：输入数据使用 [lon, lat] 坐标，输出 GeoJSON 使用 [lat, lon]。
package kriging_contour

import (
	"fmt"
	"image/color"
	"image/png"
	"math"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/duke-git/lancet/v2/maputil"
	"github.com/inkbamboo/kriging-contour/internal/contour"
	"github.com/inkbamboo/kriging-contour/internal/kriging"
	"github.com/inkbamboo/kriging-contour/internal/utils"
	"github.com/panjf2000/ants/v2"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	"github.com/paulmach/orb/planar"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"gonum.org/v1/plot/vg/vgimg"
)

// defaultLogger 包级默认日志记录器，可通过 SetLogger 替换。
var defaultLogger = utils.NewStepLogger(os.Stdout)

// SetLogger 设置全局日志记录器的输出目标。设为 nil 则抑制所有日志。
func SetLogger(writer *os.File) {
	if writer == nil {
		defaultLogger.SetLevel(utils.LogLevelError + 1) // 抑制所有级别
	} else {
		defaultLogger = utils.NewStepLogger(writer)
	}
}

// SetLogLevel 设置全局日志级别。
func SetLogLevel(level utils.LogLevel) {
	defaultLogger.SetLevel(level)
}

// GetLogger 返回当前的日志记录器，供外部调用者获取执行摘要。
func GetLogger() *utils.StepLogger {
	return defaultLogger
}

// generateGrid 基于离散数据点和边界多边形，使用克里金插值（Ordinary Kriging）生成规则网格数据。
//
// 参数:
//   - points: 原始数据点列表，X/Y 为经纬度坐标，Z 为测量值
//   - boundary: 区域边界多边形，坐标顺序为 [lat, lon]
//   - opt: 等值线参数配置，Resolution 决定网格密度
//
// 返回:
//   - grid: 插值后的网格数据，Z 矩阵行对应纬度、列对应经度
//   - err: 克里金插值失败时返回错误
//
// 插值结果会进行归一化缩放：若指定了 ContourStart/ContourEnd 则缩放到该范围，
// 否则缩放到原始数据的值域范围。
func generateGrid(points []*Point, boundary orb.Polygon, opt ContourOption) (grid *GridData, err error) {
	end := defaultLogger.BeginStep("克里金网格插值")
	defer func() {
		end(err, nil)
	}()

	minLon, maxLon := boundary.Bound().Min[1], boundary.Bound().Max[1]
	minLat, maxLat := boundary.Bound().Min[0], boundary.Bound().Max[0]
	gridX, gridY := kriging.GenerateGrid(minLon, minLat, maxLon, maxLat, opt.Resolution)

	defaultLogger.Debug("网格尺寸: %d × %d (lon: %.4f~%.4f, lat: %.4f~%.4f)",
		len(gridX), len(gridY), minLon, maxLon, minLat, maxLat)

	config := kriging.DefaultOKConfig()
	config.CoordinatesType = "euclidean"
	config.VariogramModel = "spherical"
	config.Verbose = false // 克里金内部日志关闭，用我们的统一日志

	x := make([]float64, len(points))
	y := make([]float64, len(points))
	z := make([]float64, len(points))
	for i, p := range points {
		x[i] = p.Y
		y[i] = p.X
		z[i] = p.Z
	}

	// 捕获克里金构建过程中的 panic
	var ok *kriging.OrdinaryKriging
	func() {
		defer func() {
			if r := recover(); r != nil {
				err = fmt.Errorf("克里金模型构建 panic: %v", r)
				defaultLogger.Error("克里金模型构建 panic，将使用 fallback IDW 插值: %v", r)
			}
		}()
		ok, err = kriging.NewOrdinaryKriging(x, y, z, config)
	}()

	if err != nil || ok == nil {
		err = fmt.Errorf("克里金插值失败 (将使用 IDW fallback): %w", err)
		defaultLogger.Warn("克里金插值失败，启用 IDW fallback 插值")

		// fallback: 使用逆距离加权插值 (IDW)
		grid = fallbackIDW(points, gridX, gridY)
		if grid == nil {
			return nil, fmt.Errorf("IDW fallback 插值也失败: 无法生成网格")
		}
		origMin, origMax := utils.MinMax(z...)
		gridMin, gridMax := utils.MinMax(grid.Z.RawMatrix().Data...)
		normalizeGrid(grid, gridMin, gridMax, origMin, origMax, opt)
		return grid, nil
	}

	// 捕获 ExecuteGrid 过程中的 panic
	var zvGrid *mat.Dense
	func() {
		defer func() {
			if r := recover(); r != nil {
				err = fmt.Errorf("克里金网格执行 panic: %v", r)
				defaultLogger.Error("克里金网格执行 panic，将使用 fallback IDW 插值: %v", r)
			}
		}()
		zvGrid, _ = ok.ExecuteGrid(gridX, gridY)
	}()

	if err != nil || zvGrid == nil {
		err = fmt.Errorf("克里金网格执行失败 (将使用 IDW fallback): %w", err)
		defaultLogger.Warn("克里金网格执行失败，启用 IDW fallback 插值")

		grid = fallbackIDW(points, gridX, gridY)
		if grid == nil {
			return nil, fmt.Errorf("IDW fallback 插值也失败: 无法生成网格")
		}
		origMin, origMax := utils.MinMax(z...)
		gridMin, gridMax := utils.MinMax(grid.Z.RawMatrix().Data...)
		normalizeGrid(grid, gridMin, gridMax, origMin, origMax, opt)
		return grid, nil
	}

	origMin, origMax := utils.MinMax(z...)
	gridMin, gridMax := utils.MinMax(zvGrid.RawMatrix().Data...)

	defaultLogger.Debug("Kriging 结果范围: [%.4f, %.4f], 原始范围: [%.4f, %.4f]",
		gridMin, gridMax, origMin, origMax)

	// 检查插值结果是否有效
	if math.IsNaN(gridMin) || math.IsNaN(gridMax) || math.IsInf(gridMin, 0) || math.IsInf(gridMax, 0) {
		defaultLogger.Warn("克里金插值结果包含 NaN/Inf，使用 IDW fallback")
		grid = fallbackIDW(points, gridX, gridY)
		if grid != nil {
			gMin, gMax := utils.MinMax(grid.Z.RawMatrix().Data...)
			normalizeGrid(grid, gMin, gMax, origMin, origMax, opt)
		}
		return grid, nil
	}

	grid = &GridData{
		Z: zvGrid,
		X: gridX,
		Y: gridY,
	}

	normalizeGrid(grid, gridMin, gridMax, origMin, origMax, opt)
	err = nil // 成功，清除 fallback 产生的错误
	return grid, nil
}

// normalizeGrid 对克里金插值结果进行归一化缩放。
func normalizeGrid(grid *GridData, gridMin, gridMax, origMin, origMax float64, opt ContourOption) {
	rows, cols := grid.Z.Dims()
	if gridMax-gridMin <= 1e-10 {
		return
	}
	if opt.ContourEnd-opt.ContourStart > 1e-10 {
		// 如果指定了等值线范围，缩放到等值线范围
		for i := 0; i < rows; i++ {
			for j := 0; j < cols; j++ {
				v := grid.Z.At(i, j)
				v = opt.ContourStart + (v-gridMin)/(gridMax-gridMin)*(opt.ContourEnd-opt.ContourStart)
				grid.Z.Set(i, j, v)
			}
		}
	} else if origMax-origMin > 1e-10 {
		// 否则缩放到原始数据范围
		for i := 0; i < rows; i++ {
			for j := 0; j < cols; j++ {
				v := grid.Z.At(i, j)
				v = origMin + (v-gridMin)/(gridMax-gridMin)*(origMax-origMin)
				grid.Z.Set(i, j, v)
			}
		}
	}
}

// getContourLevels 根据 ContourOption 生成等值线级别（level）列表。
// 从 ContourStart 开始，以 ContourInterval 为步长递增，不超过 ContourEnd。
// 若 ContourEnd 未指定或为零，则自动推断为数据最大值。
func getContourLevels(grid *GridData, opt ContourOption) []float64 {
	if opt.ContourInterval <= 0 {
		return nil
	}

	gridMin, gridMax := utils.MinMax(grid.Z.RawMatrix().Data...)
	start := opt.ContourStart
	end := opt.ContourEnd

	if end-start <= 1e-10 {
		// ContourEnd=0 或未指定，自动推断
		end = gridMax
	}

	var levels []float64
	for lv := start; lv <= end+1e-10; lv += opt.ContourInterval {
		if lv >= gridMin && lv <= gridMax {
			levels = append(levels, lv)
		}
	}
	return levels
}

// plotMutex 保护 gonum plot 底层 marching squares 算法中的全局缓冲区不被并发争用。
// gonum plot 的 Contour.Plot 方法内部使用了包级全局变量，并发调用会导致数据损坏。
var plotMutex sync.Mutex

// extractCoords 将 vg.Path 中的画布坐标（像素空间）逆映射为数据坐标系下的实际坐标。
// canvasSize 为画布边长（正方形），xMin/xMax/yMin/yMax 为数据域的边界范围。
func extractCoords(p vg.Path, xMin, xMax, yMin, yMax, canvasSize float64) []orb.Point {
	var pts []orb.Point
	for _, comp := range p {
		if comp.Type != vg.MoveComp && comp.Type != vg.LineComp {
			continue
		}
		x := float64(comp.Pos.X)/canvasSize*(xMax-xMin) + xMin
		y := float64(comp.Pos.Y)/canvasSize*(yMax-yMin) + yMin
		pts = append(pts, orb.Point{x, y})
	}
	return pts
}

// GenerateLines 基于离散点数据和边界多边形生成等值线（LineString）。
//
// 流程:
//  1. 调用 generateGrid 进行克里金插值生成网格数据
//  2. 使用 marching squares 算法从网格中提取等值线
//  3. 用 boundary 裁剪等值线，仅保留边界内部的部分
//  4. 合并首尾相连的线段，统一等值线方向（右手定则：高值在左）
//  5. 输出坐标转换为 [lat, lon] 格式与 boundary 对齐
//
// 参数:
//   - points: 原始数据点
//   - boundary: 区域边界多边形
//   - opt: 等值线参数（ContourStart/End/Interval/Resolution/ImagePath）
//
// 返回: GeoJSON Feature 列表，每个 Feature 为一条等值线(LineString)，
//
//	properties 含 level/is_closed/high_value_side 等字段。
//	最后一项为边界多边形 feature（properties["type"] = "boundary"）。
func GenerateLines(points []*Point, boundary orb.Polygon, opt ContourOption) []*geojson.Feature {
	// ========== 输入校验 ==========
	stats := ComputeDataStats(points)
	defaultLogger.Info("原始数据点: %d, 有效点: %d, X范围: [%.4f, %.4f], Y范围: [%.4f, %.4f], Z范围: [%.4f, %.4f], NaN: %d, Inf: %d, 重复: %d",
		stats.OriginalCount, stats.CleanedCount,
		stats.MinX, stats.MaxX, stats.MinY, stats.MaxY,
		stats.MinZ, stats.MaxZ, stats.NanCount, stats.InfCount, stats.DupCount)

	cleanPoints, warns, err := ValidatePoints(points, 3)
	for _, w := range warns {
		defaultLogger.Warn("数据校验警告: %s", w)
	}
	if err != nil {
		defaultLogger.Error("数据校验失败: %v", err)
		return nil
	}

	boundaryWarns, err := ValidateBoundary(boundary)
	for _, w := range boundaryWarns {
		defaultLogger.Warn("边界校验警告: %s", w)
	}
	if err != nil {
		defaultLogger.Error("边界校验失败: %v", err)
		return nil
	}

	optWarns, err := ValidateContourOption(opt)
	for _, w := range optWarns {
		defaultLogger.Warn("参数校验警告: %s", w)
	}
	if err != nil {
		defaultLogger.Error("参数校验失败: %v", err)
		return nil
	}

	// ========== 克里金插值 ==========
	grid, err := generateGrid(cleanPoints, boundary, opt)
	if err != nil {
		defaultLogger.Error("网格生成失败: %v", err)
		return nil
	}

	defaultLogger.Info("等值线级别数: %d", len(getContourLevels(grid, opt)))

	// ========== 等值线提取 ==========
	levelEnd := defaultLogger.BeginStep("等值线提取")
	levels := getContourLevels(grid, opt)
	if len(levels) == 0 {
		defaultLogger.Warn("没有有效的等值线级别")
		levelEnd(nil, map[string]interface{}{"warning": "无有效级别"})
		return nil
	}

	rows, cols := grid.Z.Dims()
	if rows < 2 || cols < 2 {
		defaultLogger.Error("网格尺寸过小: %d × %d", rows, cols)
		levelEnd(fmt.Errorf("网格尺寸过小"), nil)
		return nil
	}
	cg := contour.NewGrid(grid.X, grid.Y, grid.Z, cols, rows)

	const canvasSize = 2000.0

	var features []*geojson.Feature

	var wg sync.WaitGroup
	var lock sync.Mutex
	baseProperties := maputil.Merge(map[string]interface{}{
		"type": opt.ContourType,
	}, opt.Extra)
	p, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(body interface{}) {
		defer wg.Done()
		levelFeatures := generateLevelLines(cg, body.(float64), canvasSize, boundary, baseProperties)
		lock.Lock()
		features = append(features, levelFeatures...)
		defer lock.Unlock()
	})
	defer p.Release()
	// 逐层提取等值线：为每个 level 单独渲染到记录 canvas
	for _, level := range levels {
		wg.Add(1)
		_ = p.Invoke(level)
	}
	wg.Wait()

	levelEnd(nil, map[string]interface{}{"features": len(features)})
	defaultLogger.Info("等值线总数: %d", len(features))

	boundaryFeature := geojson.NewFeature(closeBoundaryPolygon(boundary))
	boundaryFeature.Properties = maputil.Merge(map[string]interface{}{"type": "boundary"}, opt.Extra)
	features = append(features, boundaryFeature)

	// 如果指定了图片路径，渲染所有等值线为一张 PNG 图片
	if opt.ImagePath != "" {
		imgEnd := defaultLogger.BeginStep("渲染等值线图片")
		renderContourImage(cg, levels, opt.ImagePath)
		imgEnd(nil, nil)
	}

	// 输出处理摘要
	defaultLogger.Summary()
	return features
}

// closeBoundaryPolygon 确保多边形所有环的首尾点一致（闭合）。
// 对于每个环，若首尾点不重合，则追加首点到末尾形成闭合环。
// 返回新的闭合多边形，不修改原始数据。
func closeBoundaryPolygon(poly orb.Polygon) orb.Polygon {
	result := make(orb.Polygon, len(poly))
	for ri, ring := range poly {
		nr := make(orb.Ring, len(ring))
		copy(nr, ring)
		if len(nr) > 0 && nr[0] != nr[len(nr)-1] {
			nr = append(nr, nr[0])
		}
		result[ri] = nr
	}
	return result
}

// generateLevelLines 为单个等值线级别生成等值线 GeoJSON Feature 列表。
//
// 流程：
//  1. 使用 gonum/plot 的 Contour.Plot 在虚拟画布上渲染该 level 的等值线
//  2. 从记录画布（Canvas）中提取 vg.Path 并逆变换为数据坐标
//  3. 用 boundary 裁剪，仅保留边界内部部分
//  4. 合并首尾相连的线段（多次迭代直到无法再合并）
//  5. 统一等值线方向：闭合环按逆时针（CCW）、开线按右手定则（高值在左）
//  6. 坐标输出前将 [lon, lat] 转换为 [lat, lon] 格式
func generateLevelLines(cg *contour.Grid, level float64, canvasSize font.Length, boundary orb.Polygon, baseProperties map[string]interface{}) []*geojson.Feature {
	var features []*geojson.Feature
	c := plotter.NewContour(cg, []float64{level}, nil)
	xMin, xMax, yMin, yMax := c.DataRange()

	rc := &contour.Canvas{}
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
	// Plot 底层 gonum marches squares 算法使用了包级全局缓冲区，并发调用会数据损坏
	plotMutex.Lock()
	c.Plot(dc, p)
	plotMutex.Unlock()

	// 收集所有坐标线
	var allCoords [][]orb.Point
	for _, path := range rc.GetPaths() {
		coords := extractCoords(path, xMin, xMax, yMin, yMax, float64(canvasSize))
		if len(coords) < 2 {
			continue
		}
		allCoords = append(allCoords, coords)
	}

	// 用 boundary 裁剪等值线，仅保留边界内部的部分
	allCoords = clipLinesToPolygon(allCoords, boundary)

	// 合并首尾相连的坐标线（误差 0.001 以内视为相连）
	tolerance := 0.01
	merged := mergeConnectedLines(allCoords, tolerance)
	for {
		oldLen := len(merged)
		merged = mergeConnectedLines(merged, tolerance)
		if len(merged) == oldLen {
			break
		}
	}

	// 计算边界周长，用于过滤过短的线（长度 < 边界周长 × 0.002）
	minLength := planar.Length(orb.LineString(boundary[0])) * 0.01
	for _, coords := range merged {
		if len(coords) < 2 {
			continue
		}
		// 过滤过短线
		if planar.Length(orb.LineString(coords)) < minLength {
			continue
		}
		// 统一处理等值线方向（右手定则: 高值永远在左侧）
		// 在坐标转换前处理（使用 [lon, lat] 坐标）
		isClosed, highSide := resolveLineOrientationCG(cg, coords, level)

		// 等值线内部计算使用 [lon, lat]，输出时转换为 [lat, lon] 与 boundary 对齐
		for i := range coords {
			coords[i][0], coords[i][1] = coords[i][1], coords[i][0]
		}

		// 坐标交换 [lon,lat] → [lat,lon] 是方向翻转变换，
		// 对闭合环：交换后 CCW → CW，需再反转回 CCW 以符合右手定则
		// 对开线：方向翻转，高值侧 left → right
		if isClosed {
			reversePoints(coords)
		} else if highSide == "left" {
			highSide = "right"
		}

		// 开线整体走向统一为逆时针（相对于边界多边形）
		// 使用边界内点判断：内点在走向左侧（cross>0）为 CCW，否则反转
		if !isClosed {
			s, e := coords[0], coords[len(coords)-1]
			rep := representativePoint(boundary)
			// 在 (x=lon, y=lat) 坐标系下计算叉积
			dx := e[1] - s[1]   // lon 方向差
			dy := e[0] - s[0]   // lat 方向差
			px := rep[1] - s[1] // 内点 lon 差
			py := rep[0] - s[0] // 内点 lat 差
			cross := dx*py - dy*px
			if cross < 0 {
				reversePoints(coords)
				if highSide == "left" {
					highSide = "right"
				} else {
					highSide = "left"
				}
			}
		}

		ls := orb.LineString(coords)
		feature := geojson.NewFeature(ls)
		feature.Properties = maputil.Merge(map[string]interface{}{}, baseProperties)
		feature.Properties["level"] = math.Round(level*1e2) / 1e2
		feature.Properties["is_closed"] = isClosed
		feature.Properties["high_value_side"] = highSide
		features = append(features, feature)
	}
	return features
}

// sampleGridCG 在 contour.Grid 中通过双线性插值采样指定坐标 (x, y) 处的值。
// 坐标系统为 [lon, lat]。若坐标超出网格范围，则取边界值。
func sampleGridCG(cg *contour.Grid, x, y float64) float64 {
	cols, rows := cg.Dims()
	if cols < 2 || rows < 2 {
		return 0
	}

	j := sort.Search(cols, func(j int) bool { return cg.X(j) >= x })
	if j >= cols {
		j = cols - 1
	}
	if j > 0 && (j >= cols || x < cg.X(j)) {
		j--
	}
	if j >= cols-1 {
		j = cols - 2
	}
	if j < 0 {
		j = 0
	}

	i := sort.Search(rows, func(i int) bool { return cg.Y(i) >= y })
	if i >= rows {
		i = rows - 1
	}
	if i > 0 && (i >= rows || y < cg.Y(i)) {
		i--
	}
	if i >= rows-1 {
		i = rows - 2
	}
	if i < 0 {
		i = 0
	}

	x0, x1 := cg.X(j), cg.X(j+1)
	y0, y1 := cg.Y(i), cg.Y(i+1)
	if x1-x0 < 1e-12 || y1-y0 < 1e-12 {
		return cg.Z(j, i)
	}

	tx := (x - x0) / (x1 - x0)
	ty := (y - y0) / (y1 - y0)
	tx = math.Max(0, math.Min(1, tx))
	ty = math.Max(0, math.Min(1, ty))

	v00 := cg.Z(j, i)
	v10 := cg.Z(j+1, i)
	v01 := cg.Z(j, i+1)
	v11 := cg.Z(j+1, i+1)

	return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
}

// resolveLineOrientationCG 调整等值线方向并确定高值侧。
//
// 对所有线条统一采用右手定则：沿线条走向，高值区始终在左侧。
// 对于闭合环，确保环为逆时针方向（CCW，有符号面积 > 0），
// 此时高值侧为 "inside"（中心值高）或 "outside"（中心值低）。
// 对于开线，通过法线方向的采样值确定方向，使高值在左，高值侧为 "left"。
//
// 注：坐标在 [lon, lat] 系统下处理。
func resolveLineOrientationCG(cg *contour.Grid, coords []orb.Point, level float64) (isClosed bool, highSide string) {
	closed := isLineClosed(coords, 0.01)
	if closed {
		n := len(coords)
		cx, cy := 0.0, 0.0
		for i := 0; i < n; i++ {
			cx += coords[i][0]
			cy += coords[i][1]
		}
		cx /= float64(n)
		cy /= float64(n)
		centerVal := sampleGridCG(cg, cx, cy)
		area := signedRingArea(coords)

		// 统一调整为逆时针方向（CCW，area > 0）
		if area < 0 {
			reversePoints(coords)
		}
		// 根据中心值与 level 的关系确定高值侧
		if centerVal > level {
			return true, "inside"
		}
		return true, "outside"
	}

	if len(coords) < 2 {
		return false, "left"
	}
	mid := len(coords) / 2
	a, b := coords[mid], coords[(mid+1)%len(coords)]
	dx := b[0] - a[0]
	dy := b[1] - a[1]
	if math.Abs(dx) < 1e-9 && math.Abs(dy) < 1e-9 {
		return false, "left"
	}

	mx := (a[0] + b[0]) / 2
	my := (a[1] + b[1]) / 2
	step := 0.001
	leftVal := sampleGridCG(cg, mx-dy*step, my+dx*step)
	rightVal := sampleGridCG(cg, mx+dy*step, my-dx*step)

	// 若右侧值大于左侧，说明当前走向高值在右，反转以保持高值在左
	if rightVal > leftVal {
		reversePoints(coords)
	}
	return false, "left"
}

// eqPoint 判断两点是否近似相等，当欧氏距离 <= tolerance 时视为相等。
func eqPoint(a, b orb.Point, tolerance float64) bool {
	return math.Hypot(a[0]-b[0], a[1]-b[1]) <= tolerance
}

// reversePoints 原地反转点序列（双指针法，O(n/2)）。
func reversePoints(pts []orb.Point) {
	for i, j := 0, len(pts)-1; i < j; i, j = i+1, j-1 {
		pts[i], pts[j] = pts[j], pts[i]
	}
}

// mergeConnectedLines 合并首尾相连的坐标线为更长的连续线段。
// 对每条未使用的线，与已有结果线逐一比较四种连接方式（首-首、尾-首、首-尾、尾-尾），
// 若匹配则合并，否则作为新线追加。tolerance 为判定两点相连的距离阈值。
func mergeConnectedLines(lines [][]orb.Point, tolerance float64) [][]orb.Point {
	if len(lines) <= 1 {
		return lines
	}
	var newLines [][]orb.Point
	newLines = append(newLines, lines[0])
	unUsedLines := lines[1:]
	for len(unUsedLines) > 0 {
		curLine := unUsedLines[0]
		unUsedLines = unUsedLines[1:]
		used := false
		for idx := 0; idx < len(newLines); idx++ {
			//首首点相同
			if eqPoint(newLines[idx][0], curLine[0], tolerance) {
				reversePoints(curLine)
				newLines[idx] = append(curLine, newLines[idx][1:]...)
				used = true
				break
			}
			//尾首点相同
			if eqPoint(newLines[idx][len(newLines[idx])-1], curLine[0], tolerance) {
				newLines[idx] = append(newLines[idx], curLine[1:]...)
				used = true
				break
			}
			//首尾点相同
			if eqPoint(newLines[idx][0], curLine[len(curLine)-1], tolerance) {
				newLines[idx] = append(curLine, newLines[idx][1:]...)
				used = true
				break
			}
			//	尾尾点相同
			if eqPoint(newLines[idx][len(newLines[idx])-1], curLine[len(curLine)-1], tolerance) {
				reversePoints(curLine)
				newLines[idx] = append(newLines[idx], curLine[1:]...)
				used = true
				break
			}
		}
		if !used {
			newLines = append(newLines, curLine)
		}
	}
	return newLines
}

// clipLinesToPolygon 用 boundary 多边形裁剪所有等值线，保留边界内部的部分。
// 注意: boundary 原始坐标顺序为 [lat, lon]，内部会统一转换为 [lon, lat] 与等值线坐标对齐。
func clipLinesToPolygon(lines [][]orb.Point, boundary orb.Polygon) [][]orb.Point {
	// boundary 原始为 [lat, lon]，转换为 [lon, lat] 与等值线坐标对齐
	normPoly := make(orb.Polygon, len(boundary))
	for ri, ring := range boundary {
		normRing := make(orb.Ring, len(ring))
		for i, p := range ring {
			normRing[i] = orb.Point{p[1], p[0]} // [lat, lon] -> [lon, lat]
		}
		normPoly[ri] = normRing
	}

	var result [][]orb.Point
	for _, line := range lines {
		clipped := clipLineToPolygon(line, normPoly)
		result = append(result, clipped...)
	}
	return result
}

// clipLineToPolygon 将单条线裁剪到 polygon 内部（均使用 [lon, lat] 坐标系）。
//
// 使用 planar.PolygonContains 判断点是否在边界内，对跨越边界的线段计算交点。
// 处理四种情况：
//   - 两点都在内部：直接添加
//   - 从内到外：找到出交点，截断当前段
//   - 从外到内：找到入交点，开始新段
//   - 两点都在外部：检查线段是否穿越多边形内部
func clipLineToPolygon(line []orb.Point, boundary orb.Polygon) [][]orb.Point {
	if len(line) < 2 {
		return nil
	}

	ring := boundary[0]
	var result [][]orb.Point
	var current []orb.Point

	prevIn := planar.PolygonContains(boundary, line[0])
	if prevIn {
		current = append(current, line[0])
	}

	for i := 1; i < len(line); i++ {
		currIn := planar.PolygonContains(boundary, line[i])

		if prevIn && currIn {
			// 两点都在内部，直接添加
			current = append(current, line[i])
		} else if prevIn && !currIn {
			// 从内部走向外部，找到出交点
			pt := lineRingFirstIntersection(line[i-1], line[i], ring)
			current = append(current, pt)
			if len(current) >= 2 {
				result = append(result, current)
			}
			current = nil
		} else if !prevIn && currIn {
			// 从外部走向内部，找到入交点
			pt := lineRingFirstIntersection(line[i-1], line[i], ring)
			current = []orb.Point{pt, line[i]}
		} else {
			// 两点都在外部，检查线段是否穿越多边形内部
			crossPts := linePolygonCrossings(line[i-1], line[i], ring)
			if len(crossPts) == 2 {
				result = append(result, []orb.Point{crossPts[0], crossPts[1]})
			}
		}

		prevIn = currIn
	}

	if len(current) >= 2 {
		result = append(result, current)
	}

	return result
}

// lineRingFirstIntersection 返回线段 a-b 与环 ring 各边的第一个交点（沿线段方向距 a 最近的交点）。
// 若未找到精确交点，则将 a 点吸附到环上的最近点作为退化处理。
func lineRingFirstIntersection(a, b orb.Point, ring orb.Ring) orb.Point {
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

	if !found {
		// 退化情况：未找到精确交点，将 a 点 snap 到环上最近点
		return snapPointToRing(a, ring)
	}
	return bestPt
}

// snapPointToRing 将点吸附（投影）到环的最近点上，返回环上的最近点坐标。
func snapPointToRing(pt orb.Point, ring orb.Ring) orb.Point {
	var best orb.Point
	bestD := math.MaxFloat64
	for i := 0; i < len(ring)-1; i++ {
		proj := closestPointOnSegment(pt, ring[i], ring[i+1])
		d := (pt[0]-proj[0])*(pt[0]-proj[0]) + (pt[1]-proj[1])*(pt[1]-proj[1])
		if d < bestD {
			bestD = d
			best = proj
		}
	}
	return best
}

// closestPointOnSegment 计算点 p 在线段 a-b 上的最近投影点。
// t 参数被钳制到 [0, 1] 区间内。
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

// linePolygonCrossings 返回线段 a-b 与环 ring 的所有交点，当 a、b 两点都在多边形外部且线段穿越内部时使用。
// 交点按距离 a 的远近排序，返回最近的两个不重复交点（构成穿入穿出对）。
func linePolygonCrossings(a, b orb.Point, ring orb.Ring) []orb.Point {
	type ptDist struct {
		pt   orb.Point
		dist float64
	}
	var crossings []ptDist

	for i := 0; i < len(ring)-1; i++ {
		pt, ok := segmentIntersection(a, b, ring[i], ring[i+1])
		if ok {
			d := (pt[0]-a[0])*(pt[0]-a[0]) + (pt[1]-a[1])*(pt[1]-a[1])
			crossings = append(crossings, ptDist{pt, d})
		}
	}

	if len(crossings) < 2 {
		return nil
	}

	// 按距 a 的距离排序，取前两个（最近的两个交点）
	sort.Slice(crossings, func(i, j int) bool {
		return crossings[i].dist < crossings[j].dist
	})

	var result []orb.Point
	// 合并距离很近的交点（容差）
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

// segmentIntersection 计算两条线段 p1-p2 和 p3-p4 的交点。
// 使用参数方程求解，t 和 u 分别表示交点在各线段上的参数位置。
// 返回 (交点, true) 如果线段相交（含端点容差），否则返回 (空点, false)。
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

// isLineClosed 判断等值线是否闭合。
// 当坐标点数 >= 3 且首尾点距离 < tol 时视为闭合。
func isLineClosed(coords []orb.Point, tol float64) bool {
	if len(coords) < 3 {
		return false
	}
	return eqPoint(coords[0], coords[len(coords)-1], tol)
}

// signedRingArea 使用鞋带公式（Shoelace formula）计算环的有符号面积。
// 正值为逆时针（CCW），负值为顺时针（CW）。
// 不含首尾重复点，通过取模实现闭合计算。
func signedRingArea(coords []orb.Point) float64 {
	if len(coords) < 3 {
		return 0
	}
	sum := 0.0
	for i := 0; i < len(coords); i++ {
		j := (i + 1) % len(coords)
		sum += coords[i][0]*coords[j][1] - coords[j][0]*coords[i][1]
	}
	return sum * 0.5
}

// renderContourImage 将克里金网格渲染为热力图，并叠加黑色等值线，输出为 PNG 图片。
// 图像尺寸固定为 1600×1200 像素，热力图使用蓝→青→绿→黄→红的渐变色带。
func renderContourImage(cg *contour.Grid, levels []float64, path string) {
	cols, rows := cg.Dims()
	if cols < 2 || rows < 2 {
		return
	}

	imgW, imgH := 1600, 1200

	xMin, xMax := cg.X(0), cg.X(cols-1)
	yMin, yMax := cg.Y(0), cg.Y(rows-1)
	if xMin > xMax {
		xMin, xMax = xMax, xMin
	}
	if yMin > yMax {
		yMin, yMax = yMax, yMin
	}

	// 生成热力图色板（256阶）
	heatColors := make([]color.Color, 256)
	for i := 0; i < 256; i++ {
		heatColors[i] = heatColorRGBA(float64(i) / 255.0)
	}

	p := plot.New()
	p.X.Min, p.X.Max = xMin, xMax
	p.Y.Min, p.Y.Max = yMin, yMax

	// 底图：热力图（所有像素着色）
	p.Add(plotter.NewHeatMap(cg, heatPalette(heatColors)))

	// 叠加：黑色等值线
	for _, lv := range levels {
		p.Add(plotter.NewContour(cg, []float64{lv}, nil))
	}

	// 渲染到同一画布
	c := vgimg.New(vg.Length(imgW), vg.Length(imgH))
	dc := draw.New(c)
	p.Draw(dc)

	f, err := os.Create(path)
	if err != nil {
		defaultLogger.Error("renderContourImage 创建文件失败: %v", err)
		return
	}
	defer f.Close()
	if err := png.Encode(f, c.Image()); err != nil {
		defaultLogger.Error("renderContourImage 编码图片失败: %v", err)
		return
	}
}

// heatPalette 简单调色板类型，实现 plot.Palette 接口的 Colors() 方法。
type heatPalette []color.Color

func (p heatPalette) Colors() []color.Color { return p }

// heatColorRGBA 返回蓝→青→绿→黄→红渐变色带的 RGBA 颜色值。
// t 取值范围 [0, 1]：
//   - 0.00-0.25: 蓝 → 青
//   - 0.25-0.50: 青 → 绿
//   - 0.50-0.75: 绿 → 黄
//   - 0.75-1.00: 黄 → 红
func heatColorRGBA(t float64) color.RGBA {
	t = math.Max(0, math.Min(1, t))
	var r, g, b uint8
	if t < 0.25 {
		s := t / 0.25
		b = 255
		g = uint8(255 * s)
	} else if t < 0.5 {
		s := (t - 0.25) / 0.25
		g = 255
		b = uint8(255 * (1 - s))
	} else if t < 0.75 {
		s := (t - 0.5) / 0.25
		g = 255
		r = uint8(255 * s)
	} else {
		s := (t - 0.75) / 0.25
		r = 255
		g = uint8(255 * (1 - s))
	}
	return color.RGBA{r, g, b, 255}
}
