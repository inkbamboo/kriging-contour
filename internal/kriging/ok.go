package kriging

import (
	"fmt"
	"math"

	"github.com/inkbamboo/kriging-contour/internal/utils"
	"gonum.org/v1/gonum/mat"
)

// OrdinaryKriging 实现普通克里金（Ordinary Kriging）插值算法。
//
// 普通克里金是一种空间插值方法，基于变差函数模型对未知点的值进行最优无偏估计。
// 该结构体保存了原始/调整后的数据坐标、变差函数模型及其参数、以及模型质量统计信息。
type OrdinaryKriging struct {
	XOrig []float64 // 原始 X 坐标
	YOrig []float64 // 原始 Y 坐标
	Z     []float64 // 观测值

	XAdjusted []float64 // 各向异性调整后的 X 坐标
	YAdjusted []float64 // 各向异性调整后的 Y 坐标

	XCenter           float64 // 坐标中心 X（用于各向异性调整）
	YCenter           float64 // 坐标中心 Y
	AnisotropyScaling float64 // 各向异性缩放因子
	AnisotropyAngle   float64 // 各向异性旋转角度（度）

	CoordinatesType string // 坐标类型："euclidean" 或 "geographic"

	VariogramModel           string        // 变差函数模型名称
	VariogramFunc            VariogramFunc // 变差函数
	VariogramModelParameters []float64     // 模型参数

	Lags         []float64 // 实验变差函数的距离滞后
	Semivariance []float64 // 实验变差函数的半方差

	Delta   []float64 // 留一交叉验证的预测误差
	Sigma   []float64 // 留一交叉验证的标准差
	Epsilon []float64 // 标准化残差 (delta/sigma)
	Q1      float64   // Q1 统计量
	Q2      float64   // Q2 统计量
	CR      float64   // cR 准则

	Verbose     bool // 是否输出详细日志
	ExactValues bool // 是否在重合点处返回精确值
	PseudoInv   bool // 是否使用伪逆而非高斯消元
}

// OKConfig 定义创建 OrdinaryKriging 实例的配置参数。
type OKConfig struct {
	VariogramModel      string         // 变差函数模型名称
	VariogramParameters interface{}    // 模型参数（[]float64 或 map[string]float64）
	VariogramFunction   VariogramFunc  // 自定义变差函数（仅 model="custom" 时有效）
	NLags               int            // 实验变差函数的分箱数
	Weight              bool           // 是否在拟合中使用加权残差
	AnisotropyScaling   float64        // 各向异性缩放因子（默认 1.0）
	AnisotropyAngle     float64        // 各向异性旋转角度（度）
	Verbose             bool           // 是否输出详细日志
	EnableStatistics    bool           // 是否计算模型质量统计
	CoordinatesType     string         // 坐标类型："euclidean" 或 "geographic"
	ExactValues         bool           // 是否在重合点处返回精确值
	PseudoInv           bool           // 是否使用伪逆
}

// DefaultOKConfig 返回默认的 OK 配置。
// 默认使用线性变差函数模型、欧几里得坐标、6 个分箱、启用精确值。
func DefaultOKConfig() OKConfig {
	return OKConfig{
		VariogramModel:    "linear",
		NLags:             6,
		AnisotropyScaling: 1.0,
		CoordinatesType:   "euclidean",
		ExactValues:       true,
	}
}

// NewOrdinaryKriging 创建并初始化一个 Ordinary Kriging 插值器。
//
// 初始化流程：
//  1. 验证输入数据的有效性
//  2. 设置变差函数模型
//  3. 应用各向异性调整
//  4. 初始化和拟合变差函数
//  5. 可选计算模型质量统计
func NewOrdinaryKriging(x, y, z []float64, cfg OKConfig) (*OrdinaryKriging, error) {
	if err := validateInputs(x, y, z, cfg); err != nil {
		return nil, err
	}

	ok := &OrdinaryKriging{
		ExactValues:     cfg.ExactValues,
		PseudoInv:       cfg.PseudoInv,
		CoordinatesType: cfg.CoordinatesType,
		Verbose:         cfg.Verbose,
	}
	ok.XOrig = append([]float64{}, x...)
	ok.YOrig = append([]float64{}, y...)
	ok.Z = append([]float64{}, z...)

	if err := ok.setupVariogramModel(cfg); err != nil {
		return nil, err
	}

	ok.setupAnisotropy(cfg)

	if ok.Verbose {
		fmt.Println("Initializing variogram model...")
	}

	vpTemp, err := MakeVariogramParameterList(cfg.VariogramModel, cfg.VariogramParameters)
	if err != nil {
		return nil, fmt.Errorf("variogram parameter error: %w", err)
	}

	ok.Lags, ok.Semivariance, ok.VariogramModelParameters, err = InitializeVariogramModel(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramModel, vpTemp, ok.VariogramFunc,
		cfg.NLags, cfg.Weight, cfg.CoordinatesType,
	)
	if err != nil {
		return nil, fmt.Errorf("variogram initialization failed: %w", err)
	}

	if ok.Verbose {
		ok.printVariogramInfo()
	}

	if cfg.EnableStatistics {
		ok.computeStatistics()
	}

	return ok, nil
}

// validateInputs 验证输入数组的长度和非空要求，以及坐标类型的有效性。
func validateInputs(x, y, z []float64, cfg OKConfig) error {
	if len(x) == 0 || len(y) == 0 || len(z) == 0 {
		return fmt.Errorf("input arrays must not be empty")
	}
	if len(x) != len(y) || len(x) != len(z) {
		return fmt.Errorf("x, y, z arrays must have the same length")
	}
	if cfg.CoordinatesType != "euclidean" && cfg.CoordinatesType != "geographic" {
		return fmt.Errorf("coordinates_type must be 'euclidean' or 'geographic'")
	}
	return nil
}

// setupVariogramModel 根据配置设置变差函数。支持内置模型和自定义函数。
func (ok *OrdinaryKriging) setupVariogramModel(cfg OKConfig) error {
	ok.VariogramModel = cfg.VariogramModel
	if fn := LookupVariogramModel(cfg.VariogramModel); fn != nil {
		ok.VariogramFunc = fn
	} else if cfg.VariogramModel == "custom" {
		if cfg.VariogramFunction == nil {
			return fmt.Errorf("must specify VariogramFunction for custom variogram model")
		}
		ok.VariogramFunc = cfg.VariogramFunction
	} else {
		return fmt.Errorf("unsupported variogram model: %s", cfg.VariogramModel)
	}
	return nil
}

// setupAnisotropy 对各向异性参数进行设置并调整坐标。
// 仅在 euclidean 坐标类型下有效；geographic 坐标忽略各向异性。
func (ok *OrdinaryKriging) setupAnisotropy(cfg OKConfig) {
	if cfg.CoordinatesType != "euclidean" {
		if cfg.AnisotropyScaling != 1.0 && ok.Verbose {
			fmt.Println("Warning: Anisotropy is not compatible with geographic coordinates. Ignoring.")
		}
		ok.AnisotropyScaling = 1.0
		ok.AnisotropyAngle = 0.0
		ok.XAdjusted = append([]float64{}, ok.XOrig...)
		ok.YAdjusted = append([]float64{}, ok.YOrig...)
		return
	}

	ok.AnisotropyScaling = cfg.AnisotropyScaling
	ok.AnisotropyAngle = cfg.AnisotropyAngle

	xMin, xMax := utils.MinMax(ok.XOrig...)
	yMin, yMax := utils.MinMax(ok.YOrig...)
	ok.XCenter = (xMax + xMin) / 2
	ok.YCenter = (yMax + yMin) / 2

	if ok.Verbose {
		fmt.Println("Adjusting data for anisotropy...")
	}

	points := make([][]float64, len(ok.XOrig))
	for i := range ok.XOrig {
		points[i] = []float64{ok.XOrig[i], ok.YOrig[i]}
	}
	adjusted := AdjustForAnisotropy(
		points,
		[]float64{ok.XCenter, ok.YCenter},
		[]float64{ok.AnisotropyScaling},
		[]float64{ok.AnisotropyAngle},
	)

	ok.XAdjusted = make([]float64, len(adjusted))
	ok.YAdjusted = make([]float64, len(adjusted))
	for i := range adjusted {
		ok.XAdjusted[i] = adjusted[i][0]
		ok.YAdjusted[i] = adjusted[i][1]
	}
}

// computeStatistics 使用留一交叉验证法计算模型质量统计（Q1、Q2、cR）。
func (ok *OrdinaryKriging) computeStatistics() {
	if ok.Verbose {
		fmt.Println("Calculating statistics on variogram model fit...")
	}
	ok.Delta, ok.Sigma, ok.Epsilon = FindStatistics(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramFunc, ok.VariogramModelParameters,
		ok.CoordinatesType, ok.PseudoInv,
	)
	ok.Q1 = CalcQ1(ok.Epsilon)
	ok.Q2 = CalcQ2(ok.Epsilon)
	ok.CR = CalcCR(ok.Q2, ok.Sigma)
	if ok.Verbose {
		fmt.Printf("Q1 = %v\nQ2 = %v\ncR = %v\n", ok.Q1, ok.Q2, ok.CR)
	}
}

// krigingMatrix 构建克里金方程组的系数矩阵 A（(n+1)×(n+1)）。
// A[i][j] = -γ(h_ij)，A[i][i] = 0，最后一行/列用于拉格朗日乘子。
func (ok *OrdinaryKriging) krigingMatrix() [][]float64 {
	n := len(ok.XAdjusted)
	n1 := n + 1

	a := make([][]float64, n1)
	for i := range a {
		a[i] = make([]float64, n1)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			dist := ok.distanceBetween(i, j)
			a[i][j] = -variogramValue(ok.VariogramFunc, ok.VariogramModelParameters, dist)
		}
		a[i][i] = 0.0
		a[n][i] = 1.0
		a[i][n] = 1.0
	}

	return a
}

// distanceBetween 计算第 i 和第 j 个已知点之间的距离。
func (ok *OrdinaryKriging) distanceBetween(i, j int) float64 {
	if ok.CoordinatesType == "euclidean" {
		return EuclideanDistance(ok.XAdjusted[i], ok.YAdjusted[i],
			ok.XAdjusted[j], ok.YAdjusted[j])
	}
	return GreatCircleDistance(ok.XAdjusted[i], ok.YAdjusted[i],
		ok.XAdjusted[j], ok.YAdjusted[j])
}

// distancesTo 计算目标点 (xpt, ypt) 到所有已知点的距离向量。
func (ok *OrdinaryKriging) distancesTo(xpt, ypt float64) []float64 {
	n := len(ok.XAdjusted)
	bd := make([]float64, n)
	if ok.CoordinatesType == "euclidean" {
		for i := 0; i < n; i++ {
			bd[i] = EuclideanDistance(ok.XAdjusted[i], ok.YAdjusted[i], xpt, ypt)
		}
	} else {
		for i := 0; i < n; i++ {
			bd[i] = GreatCircleDistance(xpt, ypt, ok.XAdjusted[i], ok.YAdjusted[i])
		}
	}
	return bd
}

// solveKrigingSystem 使用预先求逆的系数矩阵求解克里金系统。
// 返回插值值和克里金方差。
func (ok *OrdinaryKriging) solveKrigingSystem(aInv *matWrapper, bd []float64) (zvalue, sigmasq float64) {
	n := len(ok.XAdjusted)
	n1 := n + 1

	b := make([]float64, n1)
	for i := 0; i < n; i++ {
		b[i] = -variogramValue(ok.VariogramFunc, ok.VariogramModelParameters, bd[i])
		if math.Abs(bd[i]) <= eps && ok.ExactValues {
			b[i] = 0.0
		}
	}
	b[n] = 1.0

	x := aInv.MulVec(b)
	for i := 0; i < n; i++ {
		zvalue += x[i] * ok.Z[i]
	}
	for i := 0; i < n1; i++ {
		sigmasq += x[i] * (-b[i])
	}
	return
}

// executeLoop 对一系列目标点执行克里金插值。
// 先构建并求逆克里金矩阵（只需一次），然后对每个目标点求解。
// mask 为 true 的索引将被跳过（用于 masked 模式）。
func (ok *OrdinaryKriging) executeLoop(xpts, ypts []float64, mask []bool) (zvalues, sigmasq []float64) {
	npt := len(xpts)
	zvalues = make([]float64, npt)
	sigmasq = make([]float64, npt)
	a := ok.krigingMatrix()
	aInv := invertMatrix(a, ok.PseudoInv)

	for j := 0; j < npt; j++ {
		if mask != nil && mask[j] {
			continue
		}
		bd := ok.distancesTo(xpts[j], ypts[j])
		zvalues[j], sigmasq[j] = ok.solveKrigingSystem(aInv, bd)
	}
	return
}

// Execute 执行克里金插值的主方法。
//
// 支持三种模式：
//   - "grid": 在 X×Y 的笛卡尔积网格上插值
//   - "points": 在 (X[i], Y[i]) 点对序列上插值
//   - "masked": 在网格上插值但跳过 mask 为 true 的点
//
// 返回插值矩阵和方差矩阵。
func (ok *OrdinaryKriging) Execute(style string, xPoints, yPoints []float64, mask []bool) (*mat.Dense, *mat.Dense) {
	if ok.Verbose {
		fmt.Println("Executing Ordinary Kriging...")
	}

	if style != "grid" && style != "masked" && style != "points" {
		panic("style must be 'grid', 'points', or 'masked'")
	}
	xpts := append([]float64{}, xPoints...)
	ypts := append([]float64{}, yPoints...)

	var flatMask []bool
	var npt int
	nx, ny := len(xpts), len(ypts)

	switch style {
	case "grid", "masked":
		if style == "masked" {
			if len(mask) == 0 {
				panic("must specify mask when style is 'masked'")
			}
			if len(mask) != ny*nx {
				panic("mask dimensions do not match grid dimensions")
			}
			flatMask = append([]bool{}, mask...)
		}
		npt = ny * nx
		gridX := make([]float64, npt)
		gridY := make([]float64, npt)
		for i := 0; i < ny; i++ {
			for j := 0; j < nx; j++ {
				idx := i*nx + j
				gridX[idx] = xpts[j]
				gridY[idx] = ypts[i]
			}
		}
		xpts, ypts = gridX, gridY

	default:
		if nx != ny {
			panic("xPoints and yPoints must have same length for 'points' style")
		}
		npt = nx
	}
	// 对目标点同样应用各向异性调整
	if ok.CoordinatesType == "euclidean" {
		points := make([][]float64, npt)
		for i := 0; i < npt; i++ {
			points[i] = []float64{xpts[i], ypts[i]}
		}
		adjusted := AdjustForAnisotropy(points, []float64{ok.XCenter, ok.YCenter}, []float64{ok.AnisotropyScaling}, []float64{ok.AnisotropyAngle})
		for i := range adjusted {
			xpts[i], ypts[i] = adjusted[i][0], adjusted[i][1]
		}
	}

	zv, ss := ok.executeLoop(xpts, ypts, flatMask)
	if style == "grid" || style == "masked" {
		return mat.NewDense(ny, nx, zv), mat.NewDense(ny, nx, ss)
	}
	return mat.NewDense(1, npt, zv), mat.NewDense(1, npt, ss)
}

// ExecuteGrid 在 X 和 Y 坐标的笛卡尔积网格上执行克里金插值。
// 是 Execute("grid", ...) 的便捷方法。
func (ok *OrdinaryKriging) ExecuteGrid(xPoints, yPoints []float64) (*mat.Dense, *mat.Dense) {
	return ok.Execute("grid", xPoints, yPoints, nil)
}

// GetVariogramPoints 返回用于绘图或分析的拟合变差函数值。
// 返回 (距离, 拟合半方差) 对。
func (ok *OrdinaryKriging) GetVariogramPoints() ([]float64, []float64) {
	return ok.Lags, ok.VariogramFunc(ok.VariogramModelParameters, ok.Lags)
}

// GetStatistics 返回模型质量统计量 (Q1, Q2, cR)。
func (ok *OrdinaryKriging) GetStatistics() (float64, float64, float64) {
	return ok.Q1, ok.Q2, ok.CR
}

// GetEpsilonResiduals 返回标准化残差向量。
func (ok *OrdinaryKriging) GetEpsilonResiduals() []float64 {
	return ok.Epsilon
}

// SwitchVerbose 切换详细日志输出的开关状态。
func (ok *OrdinaryKriging) SwitchVerbose() {
	ok.Verbose = !ok.Verbose
}

// UpdateVariogramModel 更新变差函数模型、参数和各向异性配置。
// 更新后会自动重新拟合变差函数和计算统计量。
func (ok *OrdinaryKriging) UpdateVariogramModel(
	model string,
	params interface{},
	vfn VariogramFunc,
	nlags int, weight bool,
	anisoScaling, anisoAngle float64,
) error {

	if fn := LookupVariogramModel(model); fn != nil {
		ok.VariogramFunc = fn
	} else if model == "custom" {
		if vfn == nil {
			return fmt.Errorf("must specify variogram function for custom model")
		}
		ok.VariogramFunc = vfn
	} else {
		return fmt.Errorf("unsupported variogram model: %s", model)
	}
	ok.VariogramModel = model

	if anisoScaling != ok.AnisotropyScaling || anisoAngle != ok.AnisotropyAngle {
		if ok.CoordinatesType == "euclidean" {
			ok.AnisotropyScaling = anisoScaling
			ok.AnisotropyAngle = anisoAngle
			ok.recomputeAdjusted()
		}
	}

	vpTemp, err := MakeVariogramParameterList(model, params)
	if err != nil {
		return err
	}

	ok.Lags, ok.Semivariance, ok.VariogramModelParameters, err = InitializeVariogramModel(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramModel, vpTemp, ok.VariogramFunc,
		nlags, weight, ok.CoordinatesType,
	)
	if err != nil {
		return err
	}

	if ok.Verbose {
		ok.printVariogramInfo()
	}

	ok.computeStatistics()
	return nil
}

// recomputeAdjusted 根据当前的各向异性参数重新计算调整后的坐标。
func (ok *OrdinaryKriging) recomputeAdjusted() {
	points := make([][]float64, len(ok.XOrig))
	for i := range ok.XOrig {
		points[i] = []float64{ok.XOrig[i], ok.YOrig[i]}
	}
	adjusted := AdjustForAnisotropy(
		points,
		[]float64{ok.XCenter, ok.YCenter},
		[]float64{ok.AnisotropyScaling},
		[]float64{ok.AnisotropyAngle},
	)
	for i := range adjusted {
		ok.XAdjusted[i] = adjusted[i][0]
		ok.YAdjusted[i] = adjusted[i][1]
	}
}

// printVariogramInfo 输出当前变差函数模型的详细信息。
func (ok *OrdinaryKriging) printVariogramInfo() {
	fmt.Printf("Coordinates type: '%s'\n", ok.CoordinatesType)
	p := ok.VariogramModelParameters

	switch ok.VariogramModel {
	case "linear":
		fmt.Println("Using 'linear' Variogram Model")
		fmt.Printf("Slope: %v\n", p[0])
		fmt.Printf("Nugget: %v\n", p[1])
	case "power":
		fmt.Println("Using 'power' Variogram Model")
		fmt.Printf("Scale: %v\n", p[0])
		fmt.Printf("Exponent: %v\n", p[1])
		fmt.Printf("Nugget: %v\n", p[2])
	case "custom":
		fmt.Println("Using Custom Variogram Model")
	default:
		fmt.Printf("Using '%s' Variogram Model\n", ok.VariogramModel)
		fmt.Printf("Partial Sill: %v\n", p[0])
		fmt.Printf("Full Sill: %v\n", p[0]+p[2])
		fmt.Printf("Range: %v\n", p[1])
		fmt.Printf("Nugget: %v\n", p[2])
	}
}

// matWrapper 是一个轻量级的矩阵包装器，用于高效的矩阵-向量乘法。
// 避免完全依赖 gonum/mat 的开销。
type matWrapper struct {
	data       []float64
	rows, cols int
}

// MulVec 执行矩阵-向量乘法 y = M * b。
func (m *matWrapper) MulVec(b []float64) []float64 {
	result := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		var sum float64
		rowOffset := i * m.cols
		for j := 0; j < m.cols; j++ {
			sum += m.data[rowOffset+j] * b[j]
		}
		result[i] = sum
	}
	return result
}

// invertMatrix 对矩阵求逆。
//
// 策略：
//  1. 若 pseudoInv 为 true，使用 gonum 的伪逆
//  2. 默认使用带部分主元选取的高斯-约当消元法
//  3. 若主元过小（接近奇异），回退到伪逆
//  4. 伪逆失败时再回退到纯高斯消元法
func invertMatrix(a [][]float64, pseudoInv bool) *matWrapper {
	n := len(a)
	if pseudoInv {
		return pseudoInverse(a)
	}

	n2 := 2 * n
	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, n2)
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1.0
	}

	// 前向消元（带部分主元选取）
	for col := 0; col < n; col++ {
		pivotRow := col
		maxVal := math.Abs(aug[col][col])
		for row := col + 1; row < n; row++ {
			if v := math.Abs(aug[row][col]); v > maxVal {
				maxVal = v
				pivotRow = row
			}
		}
		if maxVal < eps {
			return pseudoInverse(a)
		}
		if pivotRow != col {
			aug[col], aug[pivotRow] = aug[pivotRow], aug[col]
		}

		pivot := aug[col][col]
		for row := col + 1; row < n; row++ {
			factor := aug[row][col] / pivot
			aug[row][col] = 0.0
			for k := col + 1; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	// 回代
	for col := n - 1; col >= 0; col-- {
		if math.Abs(aug[col][col]) < eps {
			return pseudoInverse(a)
		}
		pivot := aug[col][col]
		for k := n; k < n2; k++ {
			aug[col][k] /= pivot
		}
		for row := 0; row < col; row++ {
			factor := aug[row][col]
			for k := n; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	inv := &matWrapper{rows: n, cols: n, data: make([]float64, n*n)}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			inv.data[i*n+j] = aug[i][n+j]
		}
	}
	return inv
}

// pseudoInverse 使用 gonum 计算矩阵的伪逆（Moore-Penrose）。
// 若 gonum 求逆也失败，回退到 invertByGauss。
func pseudoInverse(a [][]float64) *matWrapper {
	n := len(a)
	if n == 0 {
		return &matWrapper{}
	}

	flat := utils.MatToFlat(a, n, n)
	matA := mat.NewDense(n, n, flat)

	var aInv mat.Dense
	if err := aInv.Inverse(matA); err != nil {
		if inv := invertByGauss(a); inv != nil {
			return inv
		}
		return &matWrapper{rows: n, cols: n, data: make([]float64, n*n)}
	}

	inv := &matWrapper{rows: n, cols: n}
	inv.data = append([]float64{}, aInv.RawMatrix().Data...)
	return inv
}

// invertByGauss 使用纯高斯-约当消元法求逆（无回退策略）。
// 作为最后的 fallback，若主元过小则返回 nil。
func invertByGauss(a [][]float64) *matWrapper {
	n := len(a)
	n2 := 2 * n

	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, n2)
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1.0
	}

	for col := 0; col < n; col++ {
		pivotRow := col
		maxVal := math.Abs(aug[col][col])
		for row := col + 1; row < n; row++ {
			if v := math.Abs(aug[row][col]); v > maxVal {
				maxVal = v
				pivotRow = row
			}
		}
		if maxVal < eps {
			return nil
		}
		if pivotRow != col {
			aug[col], aug[pivotRow] = aug[pivotRow], aug[col]
		}

		pivot := aug[col][col]
		for row := col + 1; row < n; row++ {
			factor := aug[row][col] / pivot
			aug[row][col] = 0.0
			for k := col + 1; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	for col := n - 1; col >= 0; col-- {
		pivot := aug[col][col]
		if math.Abs(pivot) < eps {
			return nil
		}
		for k := n; k < n2; k++ {
			aug[col][k] /= pivot
		}
		for row := 0; row < col; row++ {
			factor := aug[row][col]
			for k := n; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	inv := &matWrapper{rows: n, cols: n, data: make([]float64, n*n)}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			inv.data[i*n+j] = aug[i][n+j]
		}
	}
	return inv
}
