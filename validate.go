package kriging_contour

import (
	"errors"
	"fmt"
	"math"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/planar"
	"gonum.org/v1/gonum/mat"
)

// ============================================================
//  输入校验与数据清洗
// ============================================================

// ValidatePoints 校验原始数据点的合法性，自动过滤异常数据。
//
// 检查项：
//   - 点数量是否足够（至少 3 个点才能进行 kriging）
//   - 坐标是否为 NaN / Inf
//   - Z 值是否为 NaN / Inf
//   - 自动去除重复点（坐标完全相同的点）
//   - 检查点是否共线（所有点落在同一条直线上）
//
// 返回值：
//   - cleaned: 清洗后的合法数据点列表
//   - warn: 校验过程中的警告信息（非致命）
//   - err: 致命错误（数据完全不可用）
func ValidatePoints(points []*Point, minPoints int) (cleaned []*Point, warns []string, err error) {
	if minPoints <= 0 {
		minPoints = 3
	}
	if points == nil || len(points) == 0 {
		return nil, nil, errors.New("点数据为空")
	}
	if len(points) < minPoints {
		return nil, nil, fmt.Errorf("数据点数量不足 (需要至少 %d 个，当前 %d 个)", minPoints, len(points))
	}

	var clean []*Point
	seen := make(map[string]bool) // 用于去重
	nanCount, infCount, dupCount := 0, 0, 0

	for _, p := range points {
		if p == nil {
			warns = append(warns, "存在 nil 数据点，已跳过")
			continue
		}

		// 检查坐标和 Z 值的 NaN / Inf
		if math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsNaN(p.Z) {
			nanCount++
			continue
		}
		if math.IsInf(p.X, 0) || math.IsInf(p.Y, 0) || math.IsInf(p.Z, 0) {
			infCount++
			continue
		}

		// 去重（坐标完全相同的点只保留第一个）
		key := fmt.Sprintf("%.10f,%.10f", p.X, p.Y)
		if seen[key] {
			dupCount++
			continue
		}
		seen[key] = true

		clean = append(clean, &Point{X: p.X, Y: p.Y, Z: p.Z})
	}

	// 生成警告信息
	if nanCount > 0 {
		warns = append(warns, fmt.Sprintf("过滤了 %d 个 NaN 值点", nanCount))
	}
	if infCount > 0 {
		warns = append(warns, fmt.Sprintf("过滤了 %d 个 Inf 值点", infCount))
	}
	if dupCount > 0 {
		warns = append(warns, fmt.Sprintf("过滤了 %d 个重复点", dupCount))
	}

	// 检查清洗后的点数量
	if len(clean) < minPoints {
		return clean, warns, fmt.Errorf("清洗后有效数据点不足 (需要至少 %d 个，当前 %d 个)", minPoints, len(clean))
	}

	// 检查是否共线（至少 3 个点才检查）
	if len(clean) >= 3 {
		if isCollinearPoints(clean) {
			warns = append(warns, "所有数据点共线，克里金插值结果可能不准确")
		}
	}

	return clean, warns, nil
}

// isCollinearPoints 检查所有点是否近似共线。
// 使用协方差矩阵的行列式来判断：行列式接近零表示共线。
func isCollinearPoints(points []*Point) bool {
	if len(points) < 3 {
		return false
	}

	// 计算点的协方差矩阵
	var sumX, sumY float64
	for _, p := range points {
		sumX += p.X
		sumY += p.Y
	}
	meanX := sumX / float64(len(points))
	meanY := sumY / float64(len(points))

	var covXX, covYY, covXY float64
	for _, p := range points {
		dx := p.X - meanX
		dy := p.Y - meanY
		covXX += dx * dx
		covYY += dy * dy
		covXY += dx * dy
	}

	// | covXX  covXY |
	// | covXY  covYY | = covXX * covYY - covXY²
	det := covXX*covYY - covXY*covXY

	// 行列式小于阈值说明点近似共线
	// 使用相对于数据尺度的阈值
	scale := covXX + covYY
	if scale < 1e-15 {
		// 所有点的坐标基本相同，视为退化
		return true
	}

	return math.Abs(det)/scale < 1e-10
}

// ValidateBoundary 校验边界多边形的合法性。
//
// 检查项：
//   - 多边形不能为空
//   - 外环顶点数不能少于 3
//   - 外环不能自交（使用 O(n²) 线段交叉检测）
//   - 面积不能过小
func ValidateBoundary(boundary orb.Polygon) (warns []string, err error) {
	if len(boundary) == 0 {
		return nil, errors.New("边界多边形为空")
	}

	exterior := boundary[0]
	if len(exterior) < 4 {
		return nil, fmt.Errorf("边界多边形外环顶点数不足 (当前 %d 个，至少需要 3 个)", len(exterior)-1)
	}

	// 检查自交
	if isRingSelfIntersecting(exterior) {
		return nil, errors.New("边界多边形外环存在自交")
	}

	// 检查面积
	area := math.Abs(planar.Area(boundary))
	if area < 1e-8 {
		return nil, fmt.Errorf("边界多边形面积过小: %.2e", area)
	}

	return nil, nil
}

// isRingSelfIntersecting 检测多边形环是否自交。
// 使用 O(n²) 线段交叉检测，排除相邻边。
func isRingSelfIntersecting(ring orb.Ring) bool {
	n := len(ring)
	if n <= 1 {
		return false
	}
	// 去掉闭合顶点
	if ring[0] == ring[n-1] {
		n--
	}

	for i := 0; i < n; i++ {
		for j := i + 2; j < n; j++ {
			// 跳过相邻边和首尾相连的边
			if i == 0 && j == n-1 {
				continue
			}
			a, b := ring[i], ring[(i+1)%n]
			c, d := ring[j], ring[(j+1)%n]
			// 检查线段是否相交（非端点接触）
			if segmentsStrictlyIntersect(a, b, c, d) {
				return true
			}
		}
	}
	return false
}

// segmentsStrictlyIntersect 严格检查线段 AB 与 CD 是否在内部相交（排除端点接触）。
func segmentsStrictlyIntersect(a, b, c, d orb.Point) bool {
	denom := (b[0]-a[0])*(d[1]-c[1]) - (b[1]-a[1])*(d[0]-c[0])
	if math.Abs(denom) < 1e-15 {
		return false // 平行或共线
	}

	t := ((c[0]-a[0])*(d[1]-c[1]) - (c[1]-a[1])*(d[0]-c[0])) / denom
	u := -((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])) / denom

	// 严格内部相交：t 和 u 都在 (0, 1) 开区间内
	return t > 1e-12 && t < 1-1e-12 && u > 1e-12 && u < 1-1e-12
}

// ValidateContourOption 校验等值线参数配置。
func ValidateContourOption(opt ContourOption) (warns []string, err error) {
	if opt.Resolution <= 0 {
		return nil, fmt.Errorf("分辨率必须为正数，当前: %d", opt.Resolution)
	}
	if opt.Resolution < 10 {
		warns = append(warns, fmt.Sprintf("分辨率过小 (%d)，可能导致结果精度较低", opt.Resolution))
	}
	if opt.Resolution > 1000 {
		warns = append(warns, fmt.Sprintf("分辨率过大 (%d)，可能导致性能问题", opt.Resolution))
	}

	if opt.ContourInterval <= 0 {
		return nil, fmt.Errorf("等值线间隔必须为正数，当前: %f", opt.ContourInterval)
	}

	if opt.ContourStart > opt.ContourEnd && opt.ContourEnd != 0 {
		warns = append(warns, fmt.Sprintf("等值线起始值 (%f) 大于结束值 (%f)，将互换", opt.ContourStart, opt.ContourEnd))
	}

	return warns, nil
}

// ============================================================
//  统计信息收集
// ============================================================

// DataStats 数据统计信息
type DataStats struct {
	OriginalCount int
	CleanedCount  int
	NanCount      int
	InfCount      int
	DupCount      int
	MinX, MaxX    float64
	MinY, MaxY    float64
	MinZ, MaxZ    float64
	MeanZ         float64
	StdZ          float64
	IsCollinear   bool
}

// ComputeDataStats 计算数据点的统计信息。
func ComputeDataStats(points []*Point) DataStats {
	var stats DataStats
	stats.OriginalCount = len(points)

	if len(points) == 0 {
		return stats
	}

	seen := make(map[string]bool)
	var clean []*Point
	for _, p := range points {
		if p == nil {
			continue
		}
		if math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsNaN(p.Z) {
			stats.NanCount++
			continue
		}
		if math.IsInf(p.X, 0) || math.IsInf(p.Y, 0) || math.IsInf(p.Z, 0) {
			stats.InfCount++
			continue
		}
		key := fmt.Sprintf("%.10f,%.10f", p.X, p.Y)
		if seen[key] {
			stats.DupCount++
			continue
		}
		seen[key] = true
		clean = append(clean, p)
	}

	stats.CleanedCount = len(clean)
	if len(clean) == 0 {
		return stats
	}

	stats.MinX = clean[0].X
	stats.MaxX = clean[0].X
	stats.MinY = clean[0].Y
	stats.MaxY = clean[0].Y
	stats.MinZ = clean[0].Z
	stats.MaxZ = clean[0].Z
	var sumZ float64

	for _, p := range clean {
		if p.X < stats.MinX {
			stats.MinX = p.X
		}
		if p.X > stats.MaxX {
			stats.MaxX = p.X
		}
		if p.Y < stats.MinY {
			stats.MinY = p.Y
		}
		if p.Y > stats.MaxY {
			stats.MaxY = p.Y
		}
		if p.Z < stats.MinZ {
			stats.MinZ = p.Z
		}
		if p.Z > stats.MaxZ {
			stats.MaxZ = p.Z
		}
		sumZ += p.Z
	}

	stats.MeanZ = sumZ / float64(stats.CleanedCount)

	var sumSqDiff float64
	for _, p := range clean {
		diff := p.Z - stats.MeanZ
		sumSqDiff += diff * diff
	}
	stats.StdZ = math.Sqrt(sumSqDiff / float64(stats.CleanedCount))

	stats.IsCollinear = isCollinearPoints(clean)
	return stats
}

// ============================================================
//  Fallback 机制：当克里金插值失败时使用简单插值
// ============================================================

// fallbackIDW 逆距离加权插值（IDW），作为克里金插值的后备方案。
// 当克里金矩阵不可逆或拟合失败时，使用 IDW 生成网格。
func fallbackIDW(points []*Point, gridX, gridY []float64) *GridData {
	if len(points) == 0 {
		return nil
	}

	rows := len(gridY)
	cols := len(gridX)
	if rows == 0 || cols == 0 {
		return nil
	}

	// 提取坐标和值
	xs := make([]float64, len(points))
	ys := make([]float64, len(points))
	zs := make([]float64, len(points))
	for i, p := range points {
		xs[i] = p.X
		ys[i] = p.Y
		zs[i] = p.Z
	}

	// 计算数据点的空间范围用于尺度归一化
	minX, maxX := xs[0], xs[0]
	minY, maxY := ys[0], ys[0]
	for i := range xs {
		if xs[i] < minX {
			minX = xs[i]
		}
		if xs[i] > maxX {
			maxX = xs[i]
		}
		if ys[i] < minY {
			minY = ys[i]
		}
		if ys[i] > maxY {
			maxY = ys[i]
		}
	}
	scaleX := maxX - minX
	scaleY := maxY - minY
	if scaleX < 1e-10 {
		scaleX = 1.0
	}
	if scaleY < 1e-10 {
		scaleY = 1.0
	}

	// 对每个网格点进行 IDW 插值
	power := 2.0 // 距离的幂,幂越高越受近点影响
	zData := make([]float64, rows*cols)
	for r := 0; r < rows; r++ {
		for c := 0; c < cols; c++ {
			tx := gridX[c]
			ty := gridY[r]
			var weightSum, valueSum float64

			for i := range xs {
				dx := (tx - xs[i]) / scaleX
				dy := (ty - ys[i]) / scaleY
				dist := math.Sqrt(dx*dx + dy*dy)

				// 如果网格点与数据点重合，直接使用该数据点的值
				if dist < 1e-12 {
					valueSum = zs[i]
					weightSum = 1.0
					break
				}

				w := 1.0 / math.Pow(dist, power)
				weightSum += w
				valueSum += w * zs[i]
			}

			if weightSum > 0 {
				zData[r*cols+c] = valueSum / weightSum
			} else {
				// 极端情况：使用均值
				meanZ := 0.0
				for i := range zs {
					meanZ += zs[i]
				}
				zData[r*cols+c] = meanZ / float64(len(zs))
			}
		}
	}

	return &GridData{
		Z: mat.NewDense(rows, cols, zData),
		X: gridX,
		Y: gridY,
	}
}
