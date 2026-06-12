package kriging

import (
	"math"

	"github.com/inkbamboo/kriging-contour/internal/utils"
)

// GenerateGrid 根据边界多边形和分辨率创建规则插值网格的坐标轴
// 参数:
//   - boundary: 边界多边形
//   - resolution: 网格分辨率（长边方向的网格数）
//
// 返回值:
//   - []float64: X方向网格坐标列表
//   - []float64: Y方向网格坐标列表
func GenerateGrid(minX, minY, maxX, maxY float64, resolution int) ([]float64, []float64) {
	padX := (maxX - minX) * 0.02
	padY := (maxY - minY) * 0.02
	minX -= padX
	maxX += padX
	minY -= padY
	maxY += padY
	dx, dy := maxX-minX, maxY-minY
	var nx, ny int
	if dx > dy {
		nx = resolution
		ny = int(math.Max(float64(resolution)*dy/dx, 10))
	} else {
		ny = resolution
		nx = int(math.Max(float64(resolution)*dx/dy, 10))
	}

	gridX := utils.MakeLinspace(minX, maxX, nx)
	gridY := utils.MakeLinspace(minY, maxY, ny)
	return gridX, gridY
}
