package kriging

import (
	"math"

	"github.com/inkbamboo/kriging-contour/internal/utils"
)

// GenerateGrid 根据边界范围和分辨率创建规则插值网格的坐标轴。
// 在边界外增加 2% 的 padding，长边方向分辨率 = resolution，
// 短边方向按比例缩放，最小不低于 10。
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

	return utils.MakeLinspace(minX, maxX, nx), utils.MakeLinspace(minY, maxY, ny)
}
