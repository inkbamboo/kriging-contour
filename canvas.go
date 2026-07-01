// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了一个自定义的 vg.Canvas 实现，用于捕获 gonum/plot 渲染的等值线路径。
package kriging_contour

import (
	"image"
	"image/color"

	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/vg"
)

// canvas 是一个自定义的 vg.Canvas 实现，用于在 gonum/plot 渲染等值线时
// 收集 vg.Path 序列，而不进行实际的像素绘制。
// 等值线生成算法利用 gonum/plot 的 Marching Squares 实现，
// 通过此 canvas 捕获路径坐标，再转换回地理坐标。
type canvas struct {
	paths []vg.Path
}

// getPaths 返回收集到的所有 vg.Path。
func (c *canvas) getPaths() []vg.Path {
	return c.paths
}

// SetLineWidth 空实现，满足 vg.Canvas 接口。
func (c *canvas) SetLineWidth(vg.Length) {}

// SetLineDash 空实现，满足 vg.Canvas 接口。
func (c *canvas) SetLineDash([]vg.Length, vg.Length) {}

// SetColor 空实现，满足 vg.Canvas 接口。
func (c *canvas) SetColor(color.Color) {}

// Stroke 将路径追加到内部 paths 切片中，用于后续坐标提取。
func (c *canvas) Stroke(p vg.Path) { c.paths = append(c.paths, p) }

// Fill 空实现，满足 vg.Canvas 接口。
func (c *canvas) Fill(vg.Path) {}

// FillString 空实现，满足 vg.Canvas 接口。
func (c *canvas) FillString(font.Face, vg.Point, string) {}

// DrawImage 空实现，满足 vg.Canvas 接口。
func (c *canvas) DrawImage(vg.Rectangle, image.Image) {}

// Push 空实现，满足 vg.Canvas 接口。
func (c *canvas) Push() {}

// Pop 空实现，满足 vg.Canvas 接口。
func (c *canvas) Pop() {}

// Rotate 空实现，满足 vg.Canvas 接口。
func (c *canvas) Rotate(float64) {}

// Translate 空实现，满足 vg.Canvas 接口。
func (c *canvas) Translate(vg.Point) {}

// Scale 空实现，满足 vg.Canvas 接口。
func (c *canvas) Scale(float64, float64) {}
