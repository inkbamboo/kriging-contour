package contour

import (
	"gonum.org/v1/gonum/mat"
)

// Grid 将克里金插值结果 GridData 包装为 gonum plotter.GridXYZ 接口，
// 以便使用 gonum/plot 绘制等值线（Contour）和热力图（HeatMap）。
type Grid struct {
	z  *mat.Dense // 插值数据矩阵（ny 行 × nx 列）
	x  []float64  // X 轴坐标序列
	y  []float64  // Y 轴坐标序列
	nx int        // X 方向网格数
	ny int        // Y 方向网格数
}

// NewGrid 创建一个新的 Grid 实例。
func NewGrid(x, y []float64, z *mat.Dense, nx, ny int) *Grid {
	return &Grid{
		z:  z,
		x:  x,
		y:  y,
		nx: nx,
		ny: ny,
	}
}
// Dims 返回网格的维度：(列数, 行数)，满足 plotter.GridXYZ 接口。
func (g *Grid) Dims() (c, r int) { return g.nx, g.ny }

// Z 返回网格中第 c 列、第 r 行的值，满足 plotter.GridXYZ 接口。
func (g *Grid) Z(c, r int) float64 { return g.z.At(r, c) }

// X 返回第 c 列的 X 坐标，满足 plotter.GridXYZ 接口。
func (g *Grid) X(c int) float64 { return g.x[c] }

// Y 返回第 r 行的 Y 坐标，满足 plotter.GridXYZ 接口。
func (g *Grid) Y(r int) float64 { return g.y[r] }
