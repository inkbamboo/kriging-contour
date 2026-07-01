// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了地形特征提取功能，如山顶/洼地检测。
package kriging_contour

import (
	"sort"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/planar"
)

// ExtractPeaks 从插值网格中提取所有山顶的位置及高度。
//
// 识别规则：
//   - 扫描网格中所有在边界内的格点
//   - 若格点值高于其所有 8 邻域（仅比较边界内的邻居），即为局部极大值 → 山顶
//   - 边界附近的格点即使邻域不完整也算（被边界截开的半个山顶）
//   - 返回值按高度升序排列
func ExtractPeaks(grid *Grid, boundary orb.Polygon) []*Point {
	cols, rows := grid.Dims()
	if cols < 3 || rows < 3 || len(boundary) == 0 {
		return nil
	}
	// 归一化边界坐标 [Y,X] → [X,Y]
	for ri := range boundary {
		for i := range boundary[ri] {
			boundary[ri][i][0], boundary[ri][i][1] = boundary[ri][i][1], boundary[ri][i][0]
		}
	}
	// 预计算每个格点是否在边界内
	inside := make([][]bool, rows)
	for r := 0; r < rows; r++ {
		inside[r] = make([]bool, cols)
		for c := 0; c < cols; c++ {
			inside[r][c] = planar.PolygonContains(boundary, orb.Point{grid.X(c), grid.Y(r)})
		}
	}

	// 8 邻域偏移
	dirs := [][2]int{
		{-1, -1}, {-1, 0}, {-1, 1},
		{0, -1}, {0, 1},
		{1, -1}, {1, 0}, {1, 1},
	}

	var peaks []*Point
	for r := 0; r < rows; r++ {
		for c := 0; c < cols; c++ {
			if !inside[r][c] {
				continue
			}
			val := grid.Data.At(r, c)
			isPeak := true
			for _, d := range dirs {
				nr, nc := r+d[0], c+d[1]
				if nr < 0 || nr >= rows || nc < 0 || nc >= cols {
					continue
				}
				if !inside[nr][nc] {
					continue
				}
				if grid.Data.At(nr, nc) >= val {
					isPeak = false
					break
				}
			}
			if isPeak {
				peaks = append(peaks, &Point{
					X: grid.X(c),
					Y: grid.Y(r),
					Z: val,
				})
			}
		}
	}

	sort.Slice(peaks, func(i, j int) bool { return peaks[i].Z < peaks[j].Z })
	for i := range peaks {
		peaks[i].X, peaks[i].Y = peaks[i].Y, peaks[i].X
	}
	return peaks
}
