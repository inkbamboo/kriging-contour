// Package contour 提供用于等值线提取的辅助类型，包括 Grid 数据适配器和 Canvas 路径记录器。
package contour

import (
	"image"
	"image/color"

	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/vg"
)

// Canvas 是一个记录型画布，实现了 vg.Canvas 接口。
// 它仅截获 Stroke 调用以捕获等值线路径（vg.Path），其余绘制方法为空操作。
// 用于在无头模式下从 gonum/plot 的 Contour.Plot 中提取等值线数据。
type Canvas struct {
	paths []vg.Path
}

// GetPaths 返回画布记录的所有路径。
func (c *Canvas) GetPaths() []vg.Path {
	return c.paths
}
func (c *Canvas) SetLineWidth(vg.Length)                 {}
func (c *Canvas) SetLineDash([]vg.Length, vg.Length)     {}
func (c *Canvas) SetColor(color.Color)                   {}
func (c *Canvas) Stroke(p vg.Path)                       { c.paths = append(c.paths, p) }
func (c *Canvas) Fill(vg.Path)                           {}
func (c *Canvas) FillString(font.Face, vg.Point, string) {}
func (c *Canvas) DrawImage(vg.Rectangle, image.Image)    {}
func (c *Canvas) Push()                                  {}
func (c *Canvas) Pop()                                   {}
func (c *Canvas) Rotate(float64)                         {}
func (c *Canvas) Translate(vg.Point)                     {}
func (c *Canvas) Scale(float64, float64)                 {}
