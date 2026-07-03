package kriging

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// eps 数值计算中使用的零阈值。
const eps = 1.0e-10

// GreatCircleDistance 使用 Haversine 公式计算两个地理坐标点之间的大圆距离（单位：度）。
// 输入 lon1, lat1, lon2, lat2 均为角度（度），返回距离（度）。
func GreatCircleDistance(lon1, lat1, lon2, lat2 float64) float64 {
	lat1Rad := lat1 * math.Pi / 180.0
	lat2Rad := lat2 * math.Pi / 180.0
	dlon := (lon1 - lon2) * math.Pi / 180.0

	c1 := math.Cos(lat1Rad)
	s1 := math.Sin(lat1Rad)
	c2 := math.Cos(lat2Rad)
	s2 := math.Sin(lat2Rad)
	cd := math.Cos(dlon)

	return 180.0 / math.Pi * math.Atan2(
		math.Sqrt(math.Pow(c2*math.Sin(dlon), 2)+math.Pow(c1*s2-s1*c2*cd, 2)),
		s1*s2+c1*c2*cd,
	)
}

// GreatCircleDistanceVec 计算单个地理坐标点 (lon1, lat1) 到多个目标点的批量大圆距离。
func GreatCircleDistanceVec(lon1, lat1 float64, lon2, lat2 []float64) []float64 {
	n := len(lon2)
	result := make([]float64, n)
	for i := 0; i < n; i++ {
		result[i] = GreatCircleDistance(lon1, lat1, lon2[i], lat2[i])
	}
	return result
}

// EuclideanDistance 计算两点之间的欧几里得距离。
func EuclideanDistance(x1, y1, x2, y2 float64) float64 {
	dx := x1 - x2
	dy := y1 - y2
	return math.Sqrt(dx*dx + dy*dy)
}

// PairwiseEuclideanDist 计算点集中所有点对的成对欧几里得距离。
// 返回长度为 n*(n-1)/2 的上三角距离向量。
func PairwiseEuclideanDist(X, Y []float64) []float64 {
	n := len(X)
	nPairs := n * (n - 1) / 2
	dist := make([]float64, nPairs)
	idx := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dx := X[i] - X[j]
			dy := Y[i] - Y[j]
			dist[idx] = math.Sqrt(dx*dx + dy*dy)
			idx++
		}
	}
	return dist
}

// PairwiseSqEuclideanDist 计算 Z 值点集的所有成对半方差（0.5 * (zi - zj)^2）。
// 用于实验变差函数的计算。
func PairwiseSqEuclideanDist(z []float64) []float64 {
	n := len(z)
	nPairs := n * (n - 1) / 2
	sqDist := make([]float64, nPairs)
	idx := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			diff := z[i] - z[j]
			sqDist[idx] = 0.5 * diff * diff
			idx++
		}
	}
	return sqDist
}

// AdjustForAnisotropy 对点集进行各向异性调整。
// 以 center 为中心，按 angle 角度旋转，再沿 Y 方向按 scaling 缩放。
// 用于补偿数据在不同方向上的空间相关性差异。
// 支持 2D 和 3D 点集。
func AdjustForAnisotropy(X [][]float64, center, scaling, angle []float64) [][]float64 {
	nSamples := len(X)
	if nSamples == 0 {
		return X
	}
	nDim := len(X[0])

	XAdj := make([][]float64, nSamples)
	for i := 0; i < nSamples; i++ {
		XAdj[i] = make([]float64, nDim)
		for j := 0; j < nDim; j++ {
			XAdj[i][j] = X[i][j] - center[j]
		}
	}

	if nDim == 2 || nDim == 3 {
		theta := -angle[0] * math.Pi / 180.0
		cosT, sinT := math.Cos(theta), math.Sin(theta)
		s := scaling[0]
		for i := 0; i < nSamples; i++ {
			x := XAdj[i][0]*cosT - XAdj[i][1]*sinT
			y := XAdj[i][0]*sinT + XAdj[i][1]*cosT
			XAdj[i][0] = x
			XAdj[i][1] = y * s
		}
	}

	for i := 0; i < nSamples; i++ {
		for j := 0; j < nDim; j++ {
			XAdj[i][j] += center[j]
		}
	}
	return XAdj
}

// variogramParamCount 返回给定变差函数模型所需的参数数量。
func variogramParamCount(model string) int {
	switch model {
	case "linear":
		return 2 // slope, nugget
	case "power", "gaussian", "spherical", "exponential", "hole-effect":
		return 3 // psill/sill, range, nugget
	default:
		return 0
	}
}

// MakeVariogramParameterList 解析变差函数参数。
// params 可以是 []float64 或 map[string]float64 类型。
// 对于 gaussian/spherical/exponential/hole-effect 模型，会将 sill 转换为 psill（sill - nugget）。
func MakeVariogramParameterList(model string, params interface{}) ([]float64, error) {
	if params == nil {
		return nil, nil
	}

	switch v := params.(type) {
	case []float64:
		return parseListParams(model, v)
	case map[string]float64:
		return parseMapParams(model, v)
	default:
		return nil, fmt.Errorf("variogram parameters must be []float64 or map[string]float64, got %T", params)
	}
}

// parseListParams 解析 []float64 形式的变差函数参数。
// 对于 gaussian/spherical/exponential/hole-effect 模型，将 sill 转换为 psill。
func parseListParams(model string, v []float64) ([]float64, error) {
	n := variogramParamCount(model)
	if n == 0 {
		return nil, fmt.Errorf("unsupported variogram model: %s", model)
	}
	if len(v) != n {
		return nil, fmt.Errorf("%s variogram model requires exactly %d parameters, got %d", model, n, len(v))
	}
	switch model {
	case "gaussian", "spherical", "exponential", "hole-effect":
		return []float64{v[0] - v[2], v[1], v[2]}, nil
	default:
		return v, nil
	}
}

// parseMapParams 解析 map[string]float64 形式的变差函数参数。
// 支持按名称指定参数（如 "sill"、"range"、"nugget" 等）。
func parseMapParams(model string, v map[string]float64) ([]float64, error) {
	switch model {
	case "linear":
		slope, ok1 := v["slope"]
		nugget, ok2 := v["nugget"]
		if !ok1 || !ok2 {
			return nil, errors.New("linear variogram model requires 'slope' and 'nugget'")
		}
		return []float64{slope, nugget}, nil

	case "power":
		scale, ok1 := v["scale"]
		exponent, ok2 := v["exponent"]
		nugget, ok3 := v["nugget"]
		if !ok1 || !ok2 || !ok3 {
			return nil, errors.New("power variogram model requires 'scale', 'exponent', and 'nugget'")
		}
		return []float64{scale, exponent, nugget}, nil

	case "gaussian", "spherical", "exponential", "hole-effect":
		rng, ok1 := v["range"]
		nugget, ok2 := v["nugget"]
		if !ok1 || !ok2 {
			return nil, fmt.Errorf("%s variogram model requires 'range' and 'nugget'", model)
		}
		if sill, ok := v["sill"]; ok {
			return []float64{sill - nugget, rng, nugget}, nil
		}
		if psill, ok := v["psill"]; ok {
			return []float64{psill, rng, nugget}, nil
		}
		return nil, fmt.Errorf("%s variogram model requires either 'sill' or 'psill'", model)

	default:
		return nil, fmt.Errorf("unsupported variogram model: %s", model)
	}
}

// ComputeExperimentalVariogram 计算实验变差函数。
//
// 对输入数据点对的距离和半方差进行分箱统计：
//   - 使用坐标类型计算成对距离
//   - 计算 Z 值的成对半方差
//   - 将距离分 nlags 个箱
//   - 每箱输出平均距离和平均半方差
//
// 返回非空箱的 lags 和 semivariance。
func ComputeExperimentalVariogram(
	X, Y, Z []float64, nlags int, coordType string,
) (lags, semivariance []float64) {

	var d, g []float64

	if coordType == "euclidean" {
		d = PairwiseEuclideanDist(X, Y)
		g = PairwiseSqEuclideanDist(Z)
	} else {
		n := len(X)
		nPairs := n * (n - 1) / 2
		d = make([]float64, nPairs)
		g = make([]float64, nPairs)
		idx := 0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				d[idx] = GreatCircleDistance(X[i], Y[i], X[j], Y[j])
				diff := Z[i] - Z[j]
				g[idx] = 0.5 * diff * diff
				idx++
			}
		}
	}

	dmin, dmax := d[0], d[0]
	for _, v := range d {
		if v < dmin {
			dmin = v
		}
		if v > dmax {
			dmax = v
		}
	}

	binWidth := (dmax - dmin) / float64(nlags)
	bins := make([]float64, nlags+1)
	for i := 0; i < nlags; i++ {
		bins[i] = dmin + float64(i)*binWidth
	}
	bins[nlags] = dmax + 0.001

	rawLags := make([]float64, nlags)
	rawSemi := make([]float64, nlags)
	for i := 0; i < nlags; i++ {
		var sumD, sumG float64
		count := 0
		for k := range d {
			if d[k] >= bins[i] && d[k] < bins[i+1] {
				sumD += d[k]
				sumG += g[k]
				count++
			}
		}
		if count > 0 {
			rawLags[i] = sumD / float64(count)
			rawSemi[i] = sumG / float64(count)
		} else {
			rawLags[i] = math.NaN()
			rawSemi[i] = math.NaN()
		}
	}

	for i := 0; i < nlags; i++ {
		if !math.IsNaN(rawSemi[i]) {
			lags = append(lags, rawLags[i])
			semivariance = append(semivariance, rawSemi[i])
		}
	}
	return
}

// softL1 实现 Soft L1 损失函数：2 * (sqrt(1 + r^2) - 1)。
// 相比平方损失更鲁棒，对异常值不敏感。
func softL1(r float64) float64 {
	return 2.0 * (math.Sqrt(1.0+r*r) - 1.0)
}

// variogramResiduals 计算变差函数预测值与实验值的残差。
// 若 weight 为 true，使用 sigmoid 加权方案，给予中距离更高的权重。
func variogramResiduals(params, lags, semivariance []float64, vfn VariogramFunc, weight bool) []float64 {
	predicted := vfn(params, lags)
	n := len(lags)
	resid := make([]float64, n)

	if !weight {
		for i := range resid {
			resid[i] = predicted[i] - semivariance[i]
		}
		return resid
	}

	xmin, xmax := lags[0], lags[0]
	for _, v := range lags {
		if v < xmin {
			xmin = v
		} else if v > xmax {
			xmax = v
		}
	}
	dRange := xmax - xmin
	k := 2.1972 / (0.1 * dRange)
	x0 := 0.7*dRange + xmin

	weights := make([]float64, n)
	var weightSum float64
	for i := 0; i < n; i++ {
		weights[i] = 1.0 / (1.0 + math.Exp(-k*(x0-lags[i])))
		weightSum += weights[i]
	}
	for i := range weights {
		weights[i] /= weightSum
		resid[i] = (predicted[i] - semivariance[i]) * weights[i]
	}
	return resid
}

// paramBounds 存储变差函数拟合时的参数边界和初始值。
type paramBounds struct {
	lower, upper, x0 []float64 // 参数下界、上界、初始值
	semivMax         float64   // 最大半方差
	semivMin         float64   // 最小半方差
	lagsMax          float64   // 最大距离
	lagsMin          float64   // 最小距离
}

// setupParamBounds 根据模型类型和实验数据设置参数拟合的边界和初始值。
func setupParamBounds(model string, lags, semivariance []float64) paramBounds {
	nParams := variogramParamCount(model)
	b := paramBounds{
		lower: make([]float64, nParams),
		upper: make([]float64, nParams),
		x0:    make([]float64, nParams),
	}

	b.semivMax, b.semivMin = semivariance[0], semivariance[0]
	b.lagsMax, b.lagsMin = lags[0], lags[0]
	for i := range semivariance {
		if semivariance[i] > b.semivMax {
			b.semivMax = semivariance[i]
		}
		if semivariance[i] < b.semivMin {
			b.semivMin = semivariance[i]
		}
		if lags[i] > b.lagsMax {
			b.lagsMax = lags[i]
		}
		if lags[i] < b.lagsMin {
			b.lagsMin = lags[i]
		}
	}

	switch {
	case model == "linear":
		b.x0[0] = (b.semivMax - b.semivMin) / (b.lagsMax - b.lagsMin)
		b.x0[1] = b.semivMin
		b.lower[0], b.lower[1] = 0.0, 0.0
		b.upper[0], b.upper[1] = math.Inf(1), b.semivMax

	case model == "power":
		b.x0[0] = (b.semivMax - b.semivMin) / (b.lagsMax - b.lagsMin)
		b.x0[1] = 1.1
		b.x0[2] = b.semivMin
		b.lower[0], b.lower[1], b.lower[2] = 0.0, 0.001, 0.0
		b.upper[0], b.upper[1], b.upper[2] = math.Inf(1), 1.999, b.semivMax

	default:
		b.x0[0] = b.semivMax - b.semivMin
		b.x0[1] = 0.25 * b.lagsMax
		b.x0[2] = b.semivMin
		b.lower[0], b.lower[1], b.lower[2] = 0.0, 0.0, 0.0
		b.upper[0], b.upper[1], b.upper[2] = 10.0*b.semivMax, b.lagsMax, b.semivMax
	}

	return b
}

// FitVariogramModel 拟合变差函数模型参数。
//
// 使用以下优化策略：
//   - 黄金分割搜索进行一维线搜索
//   - 坐标下降法进行多维优化
//   - 多起始点策略（3 参数模型使用 5 个起始点）
//   - SoftL1 鲁棒损失函数
//   - 对过小的 range 进行矫正（避免退化为纯块金模型）
func FitVariogramModel(lags, semivariance []float64, model string, vfn VariogramFunc, weight bool) []float64 {
	nParams := variogramParamCount(model)
	bounds := setupParamBounds(model, lags, semivariance)

	clamp := func(x []float64) {
		for i := range x {
			if x[i] < bounds.lower[i] {
				x[i] = bounds.lower[i]
			}
			if !math.IsInf(bounds.upper[i], 1) && x[i] > bounds.upper[i] {
				x[i] = bounds.upper[i]
			}
		}
	}
	clamp(bounds.x0)

	objective := func(x []float64) float64 {
		r := variogramResiduals(x, lags, semivariance, vfn, weight)
		var sum float64
		for _, ri := range r {
			sum += softL1(ri)
		}
		return sum
	}

	const phi = 0.6180339887498949

	// goldenSection 对单维参数执行黄金分割搜索
	goldenSection := func(base []float64, j int, lj, uj float64) float64 {
		a, b := lj, uj
		if b-a < 1e-8 {
			return (a + b) / 2
		}
		invPhi := 1.0 - phi

		x1, x2 := a+invPhi*(b-a), a+phi*(b-a)
		p := make([]float64, nParams)
		copy(p, base)
		p[j] = x1
		f1 := objective(p)
		p[j] = x2
		f2 := objective(p)

		for b-a > 1e-6*(math.Abs(a)+math.Abs(b)) && b-a > 1e-12 {
			if f1 < f2 {
				b = x2
				x2 = x1
				f2 = f1
				x1 = a + invPhi*(b-a)
				p[j] = x1
				f1 = objective(p)
			} else {
				a = x1
				x1 = x2
				f1 = f2
				x2 = a + phi*(b-a)
				p[j] = x2
				f2 = objective(p)
			}
		}
		return (a + b) / 2
	}

	starts := [][]float64{bounds.x0}
	if nParams == 3 {
		psillUB := bounds.upper[0]
		if math.IsInf(psillUB, 1) {
			psillUB = 10.0 * (bounds.semivMax - bounds.semivMin)
		}
		rngUB := bounds.upper[1]
		if math.IsInf(rngUB, 1) {
			rngUB = bounds.lagsMax
		}
		nugUB := bounds.upper[2]
		if math.IsInf(nugUB, 1) {
			nugUB = bounds.semivMax
		}
		starts = append(starts,
			[]float64{psillUB * 0.3, rngUB * 0.15, nugUB * 0.6},
			[]float64{psillUB * 0.6, rngUB * 0.10, nugUB * 0.8},
			[]float64{psillUB * 0.1, rngUB * 0.30, nugUB * 0.9},
			[]float64{psillUB * 0.2, rngUB * 0.40, nugUB * 0.5},
		)
	}

	bestCost := math.Inf(1)
	bestResult := make([]float64, nParams)
	copy(bestResult, bounds.x0)

	// 对每个起始点执行坐标下降法
	for _, start := range starts {
		x := make([]float64, nParams)
		copy(x, start)
		clamp(x)

		for cycle := 0; cycle < 15; cycle++ {
			prevCost := objective(x)

			var order []int
			if nParams == 2 {
				order = []int{0, 1}
			} else if cycle%2 == 0 {
				order = []int{0, 1, 2}
			} else {
				order = []int{2, 1, 0}
			}

			for _, j := range order {
				lj, uj := bounds.lower[j], bounds.upper[j]
				if math.IsInf(uj, 1) {
					uj = 10.0 * (bounds.x0[0] + bounds.x0[nParams-1])
				}
				if uj <= lj || uj-lj < 1e-8 {
					continue
				}
				x[j] = goldenSection(x, j, lj, uj)
			}

			newCost := objective(x)
			if math.Abs(prevCost-newCost) < 1e-12*(1.0+math.Abs(prevCost)) {
				break
			}
		}
		clamp(x)

		if cost := objective(x); cost < bestCost {
			bestCost = cost
			copy(bestResult, x)
		}
	}

	// 对过小的 range 进行矫正（如果拟合的 range 与第一个 lag 非常接近）
	if nParams == 3 && len(lags) > 1 && bestResult[1] > 0 {
		firstLag := lags[0]
		if math.Abs(bestResult[1]-firstLag)/firstLag < 0.05 {
			corrected := firstLag * 0.75
			if corrected >= bounds.lower[1] && corrected < bestResult[1] {
				bestResult[1] = corrected
			}
		}
	}

	clamp(bestResult)
	return bestResult
}

// InitializeVariogramModel 初始化变差函数模型。
//
// 先计算实验变差函数，若提供了 modelParams 则直接使用；
// 否则调用 FitVariogramModel 自动拟合参数（custom 模型必须提供参数）。
func InitializeVariogramModel(
	X, Y, Z []float64,
	model string,
	modelParams []float64,
	vfn VariogramFunc,
	nlags int,
	weight bool,
	coordType string,
) (lags, semivariance, params []float64, err error) {
	lags, semivariance = ComputeExperimentalVariogram(X, Y, Z, nlags, coordType)

	if modelParams != nil {
		n := variogramParamCount(model)
		if n > 0 && len(modelParams) != n {
			return nil, nil, nil, fmt.Errorf("%s variogram model requires exactly %d parameters, got %d", model, n, len(modelParams))
		}
		params = modelParams
	} else {
		if model == "custom" {
			return nil, nil, nil, errors.New("must specify parameters for custom variogram model")
		}
		params = FitVariogramModel(lags, semivariance, model, vfn, weight)
	}
	return
}

// Krige 执行单点克里金插值。
//
// 构建 (n+1)×(n+1) 的克里金矩阵并求解线性方程组，返回插值值和克里金方差。
// 若目标点与已知数据点重合（距离 < eps），直接返回该点的 Z 值。
// 若矩阵求解失败，回退到平均值。
func Krige(
	X, Y, Z []float64,
	coordsX, coordsY float64,
	vfn VariogramFunc,
	vfnParams []float64,
	coordType string,
	pseudoInv bool,
) (zinterp, sigmasq float64) {

	n := len(X)
	nPlus1 := n + 1

	dMat, bd := buildKrigingDistance(X, Y, coordsX, coordsY, coordType, n)

	zeroIndex := -1
	for i, d := range bd {
		if math.Abs(d) <= eps {
			zeroIndex = i
			break
		}
	}

	aData := make([]float64, nPlus1*nPlus1)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			aData[i*nPlus1+j] = -variogramValue(vfn, vfnParams, dMat[i][j])
		}
	}
	for i := 0; i < n; i++ {
		aData[i*nPlus1+i] = 0.0
		aData[n*nPlus1+i] = 1.0
		aData[i*nPlus1+n] = 1.0
	}

	bData := make([]float64, nPlus1)
	for i := 0; i < n; i++ {
		bData[i] = -variogramValue(vfn, vfnParams, bd[i])
	}
	if zeroIndex >= 0 {
		bData[zeroIndex] = 0.0
	}
	bData[n] = 1.0

	aMat := mat.NewDense(nPlus1, nPlus1, aData)
	bVec := mat.NewVecDense(nPlus1, bData)

	var xVec mat.VecDense
	if err := xVec.SolveVec(aMat, bVec); err != nil {
		var aInv mat.Dense
		if errInv := aInv.Inverse(aMat); errInv != nil {
			var sumZ, sumW float64
			for i := 0; i < n; i++ {
				sumZ += Z[i]
				sumW += 1.0
			}
			return sumZ / sumW, 0.0
		}
		xVec.MulVec(&aInv, bVec)
	}

	for i := 0; i < n; i++ {
		zinterp += xVec.AtVec(i) * Z[i]
	}
	for i := 0; i < nPlus1; i++ {
		sigmasq += xVec.AtVec(i) * (-bData[i])
	}
	return
}

// buildKrigingDistance 构建克里金插值所需的距离矩阵。
// dMat 为已知点之间的距离矩阵，bd 为目标点到各已知点的距离向量。
func buildKrigingDistance(X, Y []float64, tx, ty float64, coordType string, n int) ([][]float64, []float64) {
	dMat := make([][]float64, n)
	bd := make([]float64, n)

	switch coordType {
	case "euclidean":
		for i := 0; i < n; i++ {
			dMat[i] = make([]float64, n)
			for j := 0; j < n; j++ {
				dMat[i][j] = EuclideanDistance(X[i], Y[i], X[j], Y[j])
			}
			bd[i] = EuclideanDistance(X[i], Y[i], tx, ty)
		}
	default:
		for i := 0; i < n; i++ {
			dMat[i] = make([]float64, n)
			for j := 0; j < n; j++ {
				dMat[i][j] = GreatCircleDistance(X[i], Y[i], X[j], Y[j])
			}
		}
		bd = GreatCircleDistanceVec(tx, ty, X, Y)
	}

	return dMat, bd
}

// FindStatistics 使用留一交叉验证法计算克里金模型的质量统计。
//
// 对每个数据点 i，用除 i 以外的所有点预测该点，计算：
//   - delta: 预测误差（Z_actual - Z_predicted）
//   - sigma: 预测标准差
//   - epsilon: 标准化残差（delta / sigma）
func FindStatistics(
	X, Y, Z []float64,
	vfn VariogramFunc,
	vfnParams []float64,
	coordType string,
	pseudoInv bool,
) (delta, sigma, epsilon []float64) {

	n := len(Z)
	allDelta := make([]float64, n)
	allSigma := make([]float64, n)

	// 预分配排除切片，避免每次循环重新分配
	xLeave := make([]float64, n-1)
	yLeave := make([]float64, n-1)
	zLeave := make([]float64, n-1)

	for i := 0; i < n; i++ {
		// 构建排除第 i 个点后的数据集
		copy(xLeave, X[:i])
		copy(xLeave[i:], X[i+1:])
		copy(yLeave, Y[:i])
		copy(yLeave[i:], Y[i+1:])
		copy(zLeave, Z[:i])
		copy(zLeave[i:], Z[i+1:])

		k, ss := Krige(
			xLeave, yLeave, zLeave,
			X[i], Y[i],
			vfn, vfnParams,
			coordType, pseudoInv,
		)
		if math.Abs(ss) < eps {
			continue
		}
		allDelta[i] = Z[i] - k
		allSigma[i] = math.Sqrt(ss)
	}

	for i := 0; i < n; i++ {
		if allSigma[i] > eps {
			delta = append(delta, allDelta[i])
			sigma = append(sigma, allSigma[i])
		}
	}

	for i := range sigma {
		epsilon = append(epsilon, delta[i]/sigma[i])
	}
	return
}

// CalcQ1 计算 Q1 统计量：标准化残差的均值绝对值。
// Q1 应接近 0，表示无系统偏差。
func CalcQ1(epsilon []float64) float64 {
	if len(epsilon) <= 1 {
		return 0
	}
	var sum float64
	for _, e := range epsilon {
		sum += e
	}
	return math.Abs(sum) / float64(len(epsilon)-1)
}

// CalcQ2 计算 Q2 统计量：标准化残差的均方。
// Q2 应接近 1，表示方差估计准确。
func CalcQ2(epsilon []float64) float64 {
	if len(epsilon) <= 1 {
		return 0
	}
	var sum float64
	for _, e := range epsilon {
		sum += e * e
	}
	return sum / float64(len(epsilon)-1)
}

// CalcCR 计算 cR 准则（组合精度准则）。
// cR = Q2 * exp(mean(log(sigma^2)))，用于综合评价模型精度。
func CalcCR(Q2 float64, sigma []float64) float64 {
	if len(sigma) == 0 {
		return 0
	}
	var sumLog float64
	for _, s := range sigma {
		sumLog += math.Log(s * s)
	}
	return Q2 * math.Exp(sumLog/float64(len(sigma)))
}
