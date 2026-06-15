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

import "math"

// VariogramFunc 是变异函数模型的函数签名。
// m 是模型参数，d 是评估距离（可能多个）。
type VariogramFunc func(m, d []float64) []float64

// variogramValue 计算单个距离值处的变异函数值，避免重复创建单元素切片。
func variogramValue(fn VariogramFunc, params []float64, dist float64) float64 {
	d := [1]float64{dist}
	return fn(params, d[:])[0]
}

// LinearVariogramModel 线性变异函数模型。
// m = [slope, nugget], γ(d) = slope * d + nugget
func LinearVariogramModel(m, d []float64) []float64 {
	slope, nugget := m[0], m[1]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = slope*dist + nugget
	}
	return result
}

// PowerVariogramModel 幂变异函数模型。
// m = [scale, exponent, nugget], γ(d) = scale * d^exponent + nugget
func PowerVariogramModel(m, d []float64) []float64 {
	scale, exponent, nugget := m[0], m[1], m[2]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = scale*math.Pow(dist, exponent) + nugget
	}
	return result
}

// GaussianVariogramModel 高斯变异函数模型。
// m = [psill, range, nugget], γ(d) = psill * (1 - exp(-d² / (range*4/7)²)) + nugget
func GaussianVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng * 4.0 / 7.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-(dist*dist)/(effRange*effRange))) + nugget
	}
	return result
}

// ExponentialVariogramModel 指数变异函数模型。
// m = [psill, range, nugget], γ(d) = psill * (1 - exp(-d / (range/3))) + nugget
func ExponentialVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// SphericalVariogramModel 球状变异函数模型。
// m = [psill, range, nugget]
// γ(d) = psill * (1.5*d/range - 0.5*(d/range)³) + nugget  (d ≤ range)
// γ(d) = psill + nugget                                    (d > range)
func SphericalVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	result := make([]float64, len(d))
	for i, dist := range d {
		if dist <= rng {
			ratio := dist / rng
			result[i] = psill*(1.5*ratio-0.5*ratio*ratio*ratio) + nugget
		} else {
			result[i] = psill + nugget
		}
	}
	return result
}

// HoleEffectVariogramModel 孔洞效应变异函数模型。
// m = [psill, range, nugget], γ(d) = psill * (1 - (1 - d/(range/3)) * exp(-d/(range/3))) + nugget
func HoleEffectVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-(1.0-dist/effRange)*math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// variogramModels 变异函数模型名称到实现的映射。
var variogramModels = map[string]VariogramFunc{
	"linear":      LinearVariogramModel,
	"power":       PowerVariogramModel,
	"gaussian":    GaussianVariogramModel,
	"spherical":   SphericalVariogramModel,
	"exponential": ExponentialVariogramModel,
	"hole-effect": HoleEffectVariogramModel,
}

// LookupVariogramModel 根据名称查找变异函数模型实现，未找到返回 nil。
func LookupVariogramModel(name string) VariogramFunc {
	return variogramModels[name]
}
