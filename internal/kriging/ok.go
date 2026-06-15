package kriging

import (
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/inkbamboo/kriging-contour/internal/utils"
	"github.com/panjf2000/ants/v2"
	"gonum.org/v1/gonum/mat"
)

// ============================================================
//  数据结构
// ============================================================

// OrdinaryKriging 提供 2D 普通克里金插值。
//
// 参考:
//   [1] P.K. Kitanidis, Introduction to Geostatistics, Cambridge University Press, 1997.
//   [2] N. Cressie, Statistics for Spatial Data, Wiley, 1993.
type OrdinaryKriging struct {
	XOrig []float64 // 原始 X 坐标
	YOrig []float64 // 原始 Y 坐标
	Z     []float64 // 观测值

	XAdjusted []float64 // 各向异性调整后的 X
	YAdjusted []float64 // 各向异性调整后的 Y

	XCenter           float64 // 中心 X
	YCenter           float64 // 中心 Y
	AnisotropyScaling float64 // 各向异性缩放
	AnisotropyAngle   float64 // 各向异性角度（度）

	CoordinatesType string // "euclidean" 或 "geographic"

	VariogramModel           string        // 模型名称
	VariogramFunc            VariogramFunc // 模型函数
	VariogramModelParameters []float64     // 拟合后参数

	Lags         []float64 // 实验变异函数 lags
	Semivariance []float64 // 实验变异函数 semivariance

	Delta   []float64 // 拟合统计: delta
	Sigma   []float64 // 拟合统计: sigma
	Epsilon []float64 // 拟合统计: epsilon
	Q1      float64   // Q1 统计量
	Q2      float64   // Q2 统计量
	CR      float64   // cR 统计量

	Verbose     bool // 是否输出进度信息
	ExactValues bool // 是否在输入位置精确插值
	PseudoInv   bool // 是否使用伪逆
}

// OKConfig 创建 OrdinaryKriging 实例的配置。
type OKConfig struct {
	VariogramModel      string        // 变异函数模型名称，默认 "linear"
	VariogramParameters interface{}   // 用户指定参数，nil 则自动拟合
	VariogramFunction   VariogramFunc // custom 模型时必须指定
	NLags               int           // 半方差图分箱数，默认 6
	Weight              bool          // 拟合时是否对近距 lag 加权
	AnisotropyScaling   float64       // 各向异性缩放，默认 1.0
	AnisotropyAngle     float64       // 各向异性角度（度，CCW），默认 0
	Verbose             bool          // 是否输出进度
	EnableStatistics    bool          // 是否计算 Q1/Q2/cR 统计量
	CoordinatesType     string        // "euclidean" 或 "geographic"，默认 "euclidean"
	ExactValues         bool          // 是否精确插值输入点（忽略 nugget），默认 true
	PseudoInv           bool          // 是否用 SVD 伪逆求解
}

// DefaultOKConfig 返回带有合理默认值的 OKConfig。
func DefaultOKConfig() OKConfig {
	return OKConfig{
		VariogramModel:    "linear",
		NLags:             6,
		AnisotropyScaling: 1.0,
		CoordinatesType:   "euclidean",
		ExactValues:       true,
	}
}

// ============================================================
//  构造函数
// ============================================================

// NewOrdinaryKriging 创建一个新的 OrdinaryKriging 实例。
func NewOrdinaryKriging(x, y, z []float64, cfg OKConfig) (*OrdinaryKriging, error) {
	if err := validateInputs(x, y, z, cfg); err != nil {
		return nil, err
	}

	ok := &OrdinaryKriging{
		XOrig:           copySlice(x),
		YOrig:           copySlice(y),
		Z:               copySlice(z),
		ExactValues:     cfg.ExactValues,
		PseudoInv:       cfg.PseudoInv,
		CoordinatesType: cfg.CoordinatesType,
		Verbose:         cfg.Verbose,
	}

	// 设置变异函数模型
	if err := ok.setupVariogramModel(cfg); err != nil {
		return nil, err
	}

	// 处理各向异性
	ok.setupAnisotropy(cfg)

	// 初始化变异函数
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

	// 计算拟合统计量
	if cfg.EnableStatistics {
		ok.computeStatistics()
	}

	return ok, nil
}

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

func (ok *OrdinaryKriging) setupAnisotropy(cfg OKConfig) {
	if cfg.CoordinatesType != "euclidean" {
		if cfg.AnisotropyScaling != 1.0 && ok.Verbose {
			fmt.Println("Warning: Anisotropy is not compatible with geographic coordinates. Ignoring.")
		}
		ok.AnisotropyScaling = 1.0
		ok.AnisotropyAngle = 0.0
		ok.XAdjusted = copySlice(ok.XOrig)
		ok.YAdjusted = copySlice(ok.YOrig)
		return
	}

	ok.AnisotropyScaling = cfg.AnisotropyScaling
	ok.AnisotropyAngle = cfg.AnisotropyAngle

	// 使用 (max+min)/2 计算中心（对齐 PyKrige）
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

// ============================================================
//  克里金矩阵构建与求解
// ============================================================

// krigingMatrix 构建 (n+1)×(n+1) 的普通克里金矩阵。
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
		a[i][i] = 0.0 // 对角线
		a[n][i] = 1.0 // Lagrange 乘数
		a[i][n] = 1.0
	}

	return a
}

// distanceBetween 返回数据集中两点 i, j 之间的距离。
func (ok *OrdinaryKriging) distanceBetween(i, j int) float64 {
	if ok.CoordinatesType == "euclidean" {
		return EuclideanDistance(ok.XAdjusted[i], ok.YAdjusted[i],
			ok.XAdjusted[j], ok.YAdjusted[j])
	}
	return GreatCircleDistance(ok.XAdjusted[i], ok.YAdjusted[i],
		ok.XAdjusted[j], ok.YAdjusted[j])
}

// distancesTo 计算数据集中所有点到目标点 (xpt, ypt) 的距离。
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

// solveKrigingSystem 使用预计算的逆矩阵求解单点的克里金系统。
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

// ============================================================
//  执行插值
// ============================================================

// executeLoop 对所有目标点逐点求解克里金系统。
func (ok *OrdinaryKriging) executeLoop(xpts, ypts []float64, mask []bool) (zvalues, sigmasq []float64) {
	npt := len(xpts)
	zvalues = make([]float64, npt)
	sigmasq = make([]float64, npt)

	a := ok.krigingMatrix()
	aInv := invertMatrix(a, ok.PseudoInv)

	var wg sync.WaitGroup
	pool, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(arg interface{}) {
		defer wg.Done()
		j := arg.(int)
		bd := ok.distancesTo(xpts[j], ypts[j])
		zvalues[j], sigmasq[j] = ok.solveKrigingSystem(aInv, bd)
	})
	defer pool.Release()

	for j := 0; j < npt; j++ {
		if mask != nil && mask[j] {
			continue
		}
		wg.Add(1)
		_ = pool.Invoke(j)
	}
	wg.Wait()
	return
}

// Execute 计算克里金插值网格及其方差。
//
// style: "grid" / "points" / "masked"
//   - grid:   xpoints 和 ypoints 定义矩形网格
//   - points: xpoints 和 ypoints 为坐标对（等长）
//   - masked: 同 grid，mask 标记跳过的点（true = 跳过）
func (ok *OrdinaryKriging) Execute(style string, xpoints, ypoints []float64, mask []bool) (*mat.Dense, *mat.Dense) {
	if ok.Verbose {
		fmt.Println("Executing Ordinary Kriging...")
	}

	if style != "grid" && style != "masked" && style != "points" {
		panic("style must be 'grid', 'points', or 'masked'")
	}

	xpts := copySlice(xpoints)
	ypts := copySlice(ypoints)

	var flatMask []bool
	var npt int
	nx, ny := len(xpts), len(ypts)

	// Meshgrid 展开
	switch style {
	case "grid", "masked":
		if style == "masked" {
			if len(mask) == 0 {
				panic("must specify mask when style is 'masked'")
			}
			if len(mask) != ny*nx {
				panic("mask dimensions do not match grid dimensions")
			}
			flatMask = copySliceBool(mask)
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

	default: // "points"
		if nx != ny {
			panic("xpoints and ypoints must have same length for 'points' style")
		}
		npt = nx
	}

	// 各向异性坐标调整
	if ok.CoordinatesType == "euclidean" {
		points := make([][]float64, npt)
		for i := 0; i < npt; i++ {
			points[i] = []float64{xpts[i], ypts[i]}
		}
		adjusted := AdjustForAnisotropy(
			points,
			[]float64{ok.XCenter, ok.YCenter},
			[]float64{ok.AnisotropyScaling},
			[]float64{ok.AnisotropyAngle},
		)
		for i := range adjusted {
			xpts[i], ypts[i] = adjusted[i][0], adjusted[i][1]
		}
	}

	// 求解 + 输出整形
	zv, ss := ok.executeLoop(xpts, ypts, flatMask)
	if style == "grid" || style == "masked" {
		return mat.NewDense(ny, nx, zv), mat.NewDense(ny, nx, ss)
	}
	return mat.NewDense(1, npt, zv), mat.NewDense(1, npt, ss)
}

// ExecuteGrid 便捷方法：以 grid 风格执行克里金插值。
func (ok *OrdinaryKriging) ExecuteGrid(xpoints, ypoints []float64) (*mat.Dense, *mat.Dense) {
	return ok.Execute("grid", xpoints, ypoints, nil)
}

// ============================================================
//  访问器
// ============================================================

// GetVariogramPoints 返回 lags 和评估后的变异函数值。
func (ok *OrdinaryKriging) GetVariogramPoints() ([]float64, []float64) {
	return ok.Lags, ok.VariogramFunc(ok.VariogramModelParameters, ok.Lags)
}

// GetStatistics 返回 Q1, Q2, cR 统计量。
func (ok *OrdinaryKriging) GetStatistics() (float64, float64, float64) {
	return ok.Q1, ok.Q2, ok.CR
}

// GetEpsilonResiduals 返回 epsilon 残差。
func (ok *OrdinaryKriging) GetEpsilonResiduals() []float64 {
	return ok.Epsilon
}

// SwitchVerbose 切换详细输出模式。
func (ok *OrdinaryKriging) SwitchVerbose() {
	ok.Verbose = !ok.Verbose
}

// ============================================================
//  更新变异函数模型
// ============================================================

// UpdateVariogramModel 更新变异函数模型和/或参数。
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

	// 更新各向异性（仅 euclidean 坐标系）
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

	// 重算统计量
	ok.computeStatistics()
	return nil
}

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

// ============================================================
//  日志输出
// ============================================================

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

// ============================================================
//  矩阵工具（内部使用）
// ============================================================

// matWrapper 是一个支持求逆的稠密矩阵包装器。
type matWrapper struct {
	data       []float64
	rows, cols int
}

// MulVec 计算矩阵与向量的乘积。
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

// invertMatrix 使用高斯消元（部分选主元）求矩阵的逆。
// 若 pseudoInv 为 true 或矩阵奇异，则使用 SVD 伪逆。
func invertMatrix(a [][]float64, pseudoInv bool) *matWrapper {
	n := len(a)
	if pseudoInv {
		return pseudoInverse(a)
	}

	// 增广矩阵 [A | I]
	n2 := 2 * n
	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, n2)
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1.0
	}

	// 前向消元（部分选主元）
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

	// 提取逆矩阵
	inv := &matWrapper{rows: n, cols: n, data: make([]float64, n*n)}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			inv.data[i*n+j] = aug[i][n+j]
		}
	}
	return inv
}

// pseudoInverse 使用 gonum 计算矩阵逆，失败时回退到高斯消元。
func pseudoInverse(a [][]float64) *matWrapper {
	n := len(a)
	if n == 0 {
		return &matWrapper{}
	}

	flat := matToFlat(a, n, n)
	matA := mat.NewDense(n, n, flat)

	var aInv mat.Dense
	if err := aInv.Inverse(matA); err != nil {
		if inv := invertByGauss(a); inv != nil {
			return inv
		}
		return &matWrapper{rows: n, cols: n}
	}

	inv := &matWrapper{rows: n, cols: n, data: make([]float64, n*n)}
	copy(inv.data, aInv.RawMatrix().Data)
	return inv
}

// invertByGauss 高斯消元（部分选主元）求逆，仅在 pseudoInverse 失败时作为最终回退。
func invertByGauss(a [][]float64) *matWrapper {
	n := len(a)
	n2 := 2 * n

	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, n2)
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1.0
	}

	// 前向消元 + 部分选主元
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
			return nil // 矩阵奇异
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

// ============================================================
//  辅助函数
// ============================================================

func copySlice[T any](s []T) []T {
	out := make([]T, len(s))
	copy(out, s)
	return out
}

func copySliceBool(s []bool) []bool {
	out := make([]bool, len(s))
	copy(out, s)
	return out
}

func matToFlat(a [][]float64, rows, cols int) []float64 {
	flat := make([]float64, rows*cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			flat[i*cols+j] = a[i][j]
		}
	}
	return flat
}
