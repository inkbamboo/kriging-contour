// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了等值线图像渲染功能，支持生成带热力图的 PNG 格式图片。
package kriging_contour

import (
	"fmt"
	"image/color"
	"image/png"
	"os"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"gonum.org/v1/plot/vg/vgimg"
)

// renderContourImage 渲染等值线图并输出为 PNG 文件。
// 首先生成热力图底色（按 level 分带着色），然后叠加等值线，
// 输出 1600×1200 像素的 PNG 图片到指定路径。
//
// 坐标系：(lat, lon)，与输出的 GeoJSON 等值线一致。
// 输入 grid 为 (lon, lat)，通过转置矩阵 + 交换坐标轴转换。
func renderContourImage(grid *Grid, levels []float64, path string) {
	cols, rows := grid.Dims()
	if cols < 2 || rows < 2 {
		return
	}

	nBands := len(levels) - 1
	if nBands <= 0 {
		return
	}

	imgW, imgH := 1600, 1200

	lonMin, lonMax := grid.X(0), grid.X(cols-1)
	latMin, latMax := grid.Y(0), grid.Y(rows-1)
	if lonMin > lonMax {
		lonMin, lonMax = lonMax, lonMin
	}
	if latMin > latMax {
		latMin, latMax = latMax, latMin
	}

	// 按 level 对网格值分带：每个格点根据其值所在的区间分配带号
	classified := mat.NewDense(rows, cols, nil)
	for r := 0; r < rows; r++ {
		for c := 0; c < cols; c++ {
			v := grid.Data.At(r, c)
			band := nBands - 1
			for i := 1; i < len(levels); i++ {
				if v <= levels[i] {
					band = i - 1
					break
				}
			}
			classified.Set(r, c, float64(band))
		}
	}

	// (lon,lat) → (lat,lon)：转置矩阵 + 交换坐标数组
	//   T 前: Data=Ny×Nx, X=lon, Y=lat
	//   T 后: Data=Nx×Ny, X=lat, Y=lon
	transClassified := &Grid{
		Data:    mat.DenseCopyOf(classified.T()),
		XCoords: grid.YCoords, // X → lat
		YCoords: grid.XCoords, // Y → lon
		Nx:      rows,
		Ny:      cols,
	}
	transGrid := &Grid{
		Data:    mat.DenseCopyOf(grid.Data.T()),
		XCoords: grid.YCoords,
		YCoords: grid.XCoords,
		Nx:      rows,
		Ny:      cols,
	}

	bandColors := levelBandColors(nBands)
	p := plot.New()
	p.X.Min, p.X.Max = latMin, latMax // X = lat
	p.Y.Min, p.Y.Max = lonMin, lonMax // Y = lon

	p.Add(plotter.NewHeatMap(transClassified, heatPalette(bandColors)))

	// 叠加每条等值线
	for _, lv := range levels {
		p.Add(plotter.NewContour(transGrid, []float64{lv}, nil))
	}

	c := vgimg.New(vg.Length(imgW), vg.Length(imgH))
	dc := draw.New(c)
	p.Draw(dc)

	f, err := os.Create(path)
	if err != nil {
		fmt.Printf("[ERROR] renderContourImage 创建文件失败: %v\n", err)
		return
	}
	defer f.Close()
	if err = png.Encode(f, c.Image()); err != nil {
		fmt.Printf("[ERROR] renderContourImage 编码图片失败: %v\n", err)
		return
	}
}

// 预定义的颜色常量，用于生成从冷色到暖色的渐变色带。
var (
	blueC   = color.RGBA{R: 0, G: 0, B: 255, A: 255}   // 蓝色（最低值）
	yellowC = color.RGBA{R: 255, G: 255, B: 0, A: 255} // 黄色
	orangeC = color.RGBA{R: 255, G: 165, B: 0, A: 255} // 橙色
	redC    = color.RGBA{R: 255, G: 0, B: 0, A: 255}   // 红色（最高值）
)

// levelBandColors 生成 n 个等值线带之间的颜色序列。
// 颜色从蓝 → 黄 → 橙 → 红渐变，用于热力图的色带映射。
// n <= 1 时返回单色（蓝色），n > 1 时生成均匀插值的渐变颜色。
func levelBandColors(n int) []color.Color {
	switch n {
	case 0:
		return nil
	case 1:
		return []color.Color{blueC}
	}

	colors := make([]color.Color, n)
	span := float64(n - 1)
	segSize := span / 3

	for i := 0; i < n; i++ {
		p := float64(i)
		switch {
		case p <= segSize:
			colors[i] = lerpColor(p/segSize, blueC, yellowC)
		case p <= 2*segSize:
			colors[i] = lerpColor((p-segSize)/segSize, yellowC, orangeC)
		default:
			colors[i] = lerpColor((p-2*segSize)/(span-2*segSize), orangeC, redC)
		}
	}
	return colors
}

// lerpColor 在两个 RGBA 颜色之间进行线性插值。
// s 为插值因子，范围 [0, 1]，0 返回 c0，1 返回 c1。
func lerpColor(s float64, c0, c1 color.RGBA) color.RGBA {
	return color.RGBA{
		R: uint8(float64(c0.R) + s*float64(c1.R-c0.R)),
		G: uint8(float64(c0.G) + s*float64(c1.G-c0.G)),
		B: uint8(float64(c0.B) + s*float64(c1.B-c0.B)),
		A: 255,
	}
}

// heatPalette 实现 plotter.Palette 接口，用于热力图的颜色映射。
type heatPalette []color.Color

// Colors 返回调色板中的颜色列表。
func (p heatPalette) Colors() []color.Color { return p }
