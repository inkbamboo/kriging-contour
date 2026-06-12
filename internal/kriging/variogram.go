// Package kriging 实现了 2D 普通克里金（Ordinary Kriging）插值算法。
//
// 本包是 PyKrige 库的 Go 语言移植，支持多种变异函数模型：
//   - linear: 线性模型
//   - power: 幂模型
//   - gaussian: 高斯模型
//   - spherical: 球状模型（默认）
//   - exponential: 指数模型
//   - hole-effect: 孔洞效应模型
//
// 支持欧氏距离和地理距离（大圆距离）两种坐标系。
package kriging

import (
	"math"
)

// VariogramFunc is the function signature for a variogram model.
// m contains the model parameters, d is the distance(s) at which to evaluate.
type VariogramFunc func(m []float64, d []float64) []float64

// LinearVariogramModel implements the linear variogram model.
// m = [slope, nugget]
// γ(d) = slope * d + nugget
func LinearVariogramModel(m []float64, d []float64) []float64 {
	slope := m[0]
	nugget := m[1]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = slope*dist + nugget
	}
	return result
}

// PowerVariogramModel implements the power variogram model.
// m = [scale, exponent, nugget]
// γ(d) = scale * d^exponent + nugget
func PowerVariogramModel(m []float64, d []float64) []float64 {
	scale := m[0]
	exponent := m[1]
	nugget := m[2]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = scale*math.Pow(dist, exponent) + nugget
	}
	return result
}

// GaussianVariogramModel implements the gaussian variogram model.
// m = [psill, range, nugget]
// γ(d) = psill * (1 - exp(-d² / (range*4/7)²)) + nugget
func GaussianVariogramModel(m []float64, d []float64) []float64 {
	psill := m[0]
	rangeVal := m[1]
	nugget := m[2]
	// effective range = range * 4/7
	effRange := rangeVal * 4.0 / 7.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-(dist*dist)/(effRange*effRange))) + nugget
	}
	return result
}

// ExponentialVariogramModel implements the exponential variogram model.
// m = [psill, range, nugget]
// γ(d) = psill * (1 - exp(-d / (range/3))) + nugget
func ExponentialVariogramModel(m []float64, d []float64) []float64 {
	psill := m[0]
	rangeVal := m[1]
	nugget := m[2]
	// effective range = range / 3
	effRange := rangeVal / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// SphericalVariogramModel implements the spherical variogram model.
// m = [psill, range, nugget]
// γ(d) = psill * (1.5*d/range - 0.5*(d/range)³) + nugget  for d <= range
// γ(d) = psill + nugget                                    for d > range
func SphericalVariogramModel(m []float64, d []float64) []float64 {
	psill := m[0]
	rangeVal := m[1]
	nugget := m[2]
	result := make([]float64, len(d))
	for i, dist := range d {
		if dist <= rangeVal {
			ratio := dist / rangeVal
			result[i] = psill*(1.5*ratio-0.5*ratio*ratio*ratio) + nugget
		} else {
			result[i] = psill + nugget
		}
	}
	return result
}

// HoleEffectVariogramModel implements the hole-effect variogram model.
// m = [psill, range, nugget]
// γ(d) = psill * (1 - (1 - d/(range/3)) * exp(-d/(range/3))) + nugget
func HoleEffectVariogramModel(m []float64, d []float64) []float64 {
	psill := m[0]
	rangeVal := m[1]
	nugget := m[2]
	effRange := rangeVal / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-(1.0-dist/effRange)*math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// VariogramModelMap maps variogram model names to their implementations.
var VariogramModelMap = map[string]VariogramFunc{
	"linear":      LinearVariogramModel,
	"power":       PowerVariogramModel,
	"gaussian":    GaussianVariogramModel,
	"spherical":   SphericalVariogramModel,
	"exponential": ExponentialVariogramModel,
	"hole-effect": HoleEffectVariogramModel,
}
