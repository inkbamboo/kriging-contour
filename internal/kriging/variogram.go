package kriging

import "math"

// VariogramFunc 定义变差函数类型。
// 参数 m 为模型参数数组，d 为距离数组，返回对应距离下的半方差值。
type VariogramFunc func(m, d []float64) []float64

// variogramValue 计算单个距离点的变差函数值。
// 是 VariogramFunc 的单值版本辅助函数。
func variogramValue(fn VariogramFunc, params []float64, dist float64) float64 {
	d := [1]float64{dist}
	return fn(params, d[:])[0]
}

// LinearVariogramModel 线性变差函数模型。
// 参数: m[0]=slope（斜率）, m[1]=nugget（块金值）。
// 公式: γ(d) = slope * d + nugget
func LinearVariogramModel(m, d []float64) []float64 {
	slope, nugget := m[0], m[1]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = slope*dist + nugget
	}
	return result
}

// PowerVariogramModel 幂变差函数模型。
// 参数: m[0]=scale（比例）, m[1]=exponent（指数，范围 0<exp<2）, m[2]=nugget（块金值）。
// 公式: γ(d) = scale * d^exponent + nugget
func PowerVariogramModel(m, d []float64) []float64 {
	scale, exponent, nugget := m[0], m[1], m[2]
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = scale*math.Pow(dist, exponent) + nugget
	}
	return result
}

// GaussianVariogramModel 高斯变差函数模型。
// 参数: m[0]=psill（偏基台值）, m[1]=range（范围）, m[2]=nugget（块金值）。
// 有效范围 effRange = range * 4/7（约 57% 的范围处达到 95% 的基台值）。
// 公式: γ(d) = psill * (1 - exp(-d²/effRange²)) + nugget
func GaussianVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng * 4.0 / 7.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-(dist*dist)/(effRange*effRange))) + nugget
	}
	return result
}

// ExponentialVariogramModel 指数变差函数模型。
// 参数: m[0]=psill（偏基台值）, m[1]=range（范围）, m[2]=nugget（块金值）。
// 有效范围 effRange = range / 3（在 range 处达到约 95% 的基台值）。
// 公式: γ(d) = psill * (1 - exp(-d/effRange)) + nugget
func ExponentialVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// SphericalVariogramModel 球状变差函数模型（最常用的地统计模型之一）。
// 参数: m[0]=psill（偏基台值）, m[1]=range（范围）, m[2]=nugget（块金值）。
// 公式:
//
//	若 d <= range: γ(d) = psill * (1.5*r - 0.5*r³) + nugget  (r = d/range)
//	若 d >  range: γ(d) = psill + nugget
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

// HoleEffectVariogramModel 孔穴效应变差函数模型。
// 适用于存在周期性波动（孔穴效应）的数据。
// 参数: m[0]=psill（偏基台值）, m[1]=range（范围）, m[2]=nugget（块金值）。
// 有效范围 effRange = range / 3。
// 公式: γ(d) = psill * (1 - (1 - d/effRange) * exp(-d/effRange)) + nugget
func HoleEffectVariogramModel(m, d []float64) []float64 {
	psill, rng, nugget := m[0], m[1], m[2]
	effRange := rng / 3.0
	result := make([]float64, len(d))
	for i, dist := range d {
		result[i] = psill*(1.0-(1.0-dist/effRange)*math.Exp(-dist/effRange)) + nugget
	}
	return result
}

// variogramModels 内置变差函数模型名称到函数的映射。
var variogramModels = map[string]VariogramFunc{
	"linear":      LinearVariogramModel,
	"power":       PowerVariogramModel,
	"gaussian":    GaussianVariogramModel,
	"spherical":   SphericalVariogramModel,
	"exponential": ExponentialVariogramModel,
	"hole-effect": HoleEffectVariogramModel,
}

// LookupVariogramModel 根据名称查找内置变差函数模型。
// 返回 nil 表示模型不存在（需使用自定义函数）。
func LookupVariogramModel(name string) VariogramFunc {
	return variogramModels[name]
}
