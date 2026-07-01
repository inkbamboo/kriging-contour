// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了数据校验、参数验证和统计信息计算功能。
package kriging_contour

import (
	"errors"
	"fmt"
	"math"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/planar"
)

// ValidatePoints 对输入数据点进行清洗和校验。
// 清洗过程包括：过滤 nil 点、NaN 值、Inf 值，按 (X, Y) 坐标去重。
// 同时检测点是否共线（共线数据不适合克里金插值）。
// minPoints 为最小有效点数要求，若 <= 0 则默认为 3。
// 返回清洗后的有效点切片和错误信息。
func ValidatePoints(points []*Point, minPoints int) (cleaned []*Point, err error) {
	if minPoints <= 0 {
		minPoints = 3
	}
	if points == nil || len(points) == 0 {
		return nil, errors.New("点数据为空")
	}
	if len(points) < minPoints {
		return nil, fmt.Errorf("数据点数量不足 (需要至少 %d 个，当前 %d 个)", minPoints, len(points))
	}

	var clean []*Point
	seen := make(map[string]bool)
	nanCount, infCount, dupCount := 0, 0, 0

	for _, p := range points {
		if p == nil {
			fmt.Println("[WARN] 数据校验警告: 存在 nil 数据点，已跳过")
			continue
		}

		if math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsNaN(p.Z) {
			nanCount++
			continue
		}
		if math.IsInf(p.X, 0) || math.IsInf(p.Y, 0) || math.IsInf(p.Z, 0) {
			infCount++
			continue
		}

		key := fmt.Sprintf("%.10f,%.10f", p.X, p.Y)
		if seen[key] {
			dupCount++
			continue
		}
		seen[key] = true

		clean = append(clean, &Point{X: p.X, Y: p.Y, Z: p.Z})
	}

	if nanCount > 0 {
		fmt.Printf("[WARN] 数据校验警告: 过滤了 %d 个 NaN 值点\n", nanCount)
	}
	if infCount > 0 {
		fmt.Printf("[WARN] 数据校验警告: 过滤了 %d 个 Inf 值点\n", infCount)
	}
	if dupCount > 0 {
		fmt.Printf("[WARN] 数据校验警告: 过滤了 %d 个重复点\n", dupCount)
	}

	if len(clean) < minPoints {
		return clean, fmt.Errorf("清洗后有效数据点不足 (需要至少 %d 个，当前 %d 个)", minPoints, len(clean))
	}

	if len(clean) >= 3 {
		if isCollinearPoints(clean) {
			fmt.Println("[WARN] 数据校验警告: 所有数据点共线，克里金插值结果可能不准确")
		}
	}

	return clean, nil
}

// isCollinearPoints 通过 PCA 思想检测所有数据点是否近似共线。
// 计算 X、Y 坐标的协方差矩阵行列式，若行列式相对缩放因子极小，
// 则判断为共线。共线的点无法进行有意义的二维克里金插值。
func isCollinearPoints(points []*Point) bool {
	if len(points) < 3 {
		return false
	}

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

	det := covXX*covYY - covXY*covXY
	scale := covXX + covYY
	if scale < 1e-15 {
		return true
	}

	return math.Abs(det)/scale < 1e-10
}

// ValidateBoundary 校验边界多边形的有效性。
// 检查内容包括：多边形是否为空、外环顶点数是否足够、外环是否自交、
// 面积是否过小。边界多边形用于裁剪等值线和限制插值范围。
func ValidateBoundary(boundary orb.Polygon) error {
	if len(boundary) == 0 {
		return errors.New("边界多边形为空")
	}

	exterior := boundary[0]
	if len(exterior) < 4 {
		return fmt.Errorf("边界多边形外环顶点数不足 (当前 %d 个，至少需要 3 个)", len(exterior)-1)
	}

	if isRingSelfIntersecting(exterior) {
		return errors.New("边界多边形外环存在自交")
	}

	area := math.Abs(planar.Area(boundary))
	if area < 1e-8 {
		return fmt.Errorf("边界多边形面积过小: %.2e", area)
	}

	return nil
}

// isRingSelfIntersecting 检测多边形环是否存在自交。
// 遍历所有非相邻线段对，使用严格相交判断算法检测。
// 相邻边和通过首尾相连的边对不视为自交。
func isRingSelfIntersecting(ring orb.Ring) bool {
	n := len(ring)
	if n <= 1 {
		return false
	}
	if ring[0] == ring[n-1] {
		n--
	}

	for i := 0; i < n; i++ {
		for j := i + 2; j < n; j++ {
			if i == 0 && j == n-1 {
				continue
			}
			a, b := ring[i], ring[(i+1)%n]
			c, d := ring[j], ring[(j+1)%n]
			if segmentsStrictlyIntersect(a, b, c, d) {
				return true
			}
		}
	}
	return false
}

// segmentsStrictlyIntersect 判断两条线段是否严格相交（不包括端点重合的情况）。
// 使用参数方法计算交点参数 t 和 u，仅当二者均在开区间 (0, 1) 时返回 true。
func segmentsStrictlyIntersect(a, b, c, d orb.Point) bool {
	denom := (b[0]-a[0])*(d[1]-c[1]) - (b[1]-a[1])*(d[0]-c[0])
	if math.Abs(denom) < 1e-15 {
		return false
	}

	t := ((c[0]-a[0])*(d[1]-c[1]) - (c[1]-a[1])*(d[0]-c[0])) / denom
	u := -((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])) / denom

	return t > 1e-12 && t < 1-1e-12 && u > 1e-12 && u < 1-1e-12
}

// ValidateContourOption 校验等值线配置参数。
// 校验 Resolution 范围，校验 LevelList 或 ContourInterval 必须指定一项。
// LevelList 的实际计算由 GetLevelList 完成，调用方负责在合适的时机补全。
func ValidateContourOption(opt *ContourOption) error {
	if opt.Resolution <= 0 {
		return fmt.Errorf("分辨率必须为正数，当前: %d", opt.Resolution)
	}
	if opt.Resolution < 10 {
		fmt.Printf("[WARN] 参数校验警告: 分辨率过小 (%d)，可能导致结果精度较低\n", opt.Resolution)
	}
	if opt.Resolution > 1000 {
		fmt.Printf("[WARN] 参数校验警告: 分辨率过大 (%d)，可能导致性能问题\n", opt.Resolution)
	}

	if len(opt.LevelList) == 0 && opt.ContourInterval <= 0 {
		return fmt.Errorf("必须指定 level_list 或 contour_start/contour_end/contour_interval")
	}

	return nil
}

// DataStats 存储输入数据的统计信息。
// 包含数据清洗前后的计数、各坐标维度的范围、Z 值的均值/标准差，
// 以及是否共线的判断结果。
type DataStats struct {
	OriginalCount int     // 原始数据点数
	CleanedCount  int     // 清洗后有效数据点数
	NanCount      int     // NaN 值数量
	InfCount      int     // Inf 值数量
	DupCount      int     // 重复点数量
	MinX, MaxX    float64 // X 坐标范围
	MinY, MaxY    float64 // Y 坐标范围
	MinZ, MaxZ    float64 // Z 值范围
	MeanZ         float64 // Z 值均值
	StdZ          float64 // Z 值标准差
	IsCollinear   bool    // 数据是否近似共线
}

// ComputeDataStats 计算输入数据点的统计信息（不修改原始数据）。
// 执行与 ValidatePoints 相同的 NaN/Inf/重复点过滤，但不返回清洗后的数据，
// 只返回统计数据。用于在 GenerateLines 中输出数据概况。
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
