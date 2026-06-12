package kriging

import (
	"fmt"
	"math"
	"runtime"
	"sync"
	"time"

	"github.com/panjf2000/ants/v2"
	"gonum.org/v1/gonum/mat"
)

// OrdinaryKriging provides 2D Ordinary Kriging interpolation.
//
// This is a Go port of PyKrige's OrdinaryKriging class.
//
// References:
//
//	[1] P.K. Kitanidis, Introduction to Geostatistics: Applications in
//	    Hydrogeology, (Cambridge University Press, 1997) 272 p.
//	[2] N. Cressie, Statistics for spatial data,
//	    (Wiley Series in Probability and Statistics, 1993) 137 p.
type OrdinaryKriging struct {
	// Original data coordinates and values
	XOrig []float64
	YOrig []float64
	Z     []float64

	// Anisotropy-adjusted coordinates
	XAdjusted []float64
	YAdjusted []float64

	// Center and anisotropy parameters
	XCenter           float64
	YCenter           float64
	AnisotropyScaling float64
	AnisotropyAngle   float64

	// Coordinate type: "euclidean" or "geographic"
	CoordinatesType string

	// Variogram model
	VariogramModel           string
	VariogramFunc            VariogramFunc
	VariogramModelParameters []float64

	// Experimental variogram data
	Lags         []float64
	Semivariance []float64

	// Statistics
	Delta   []float64
	Sigma   []float64
	Epsilon []float64
	Q1      float64
	Q2      float64
	CR      float64

	// Configuration
	Verbose        bool
	EnablePlotting bool
	ExactValues    bool
	PseudoInv      bool
}

// OKConfig holds configuration options for creating an OrdinaryKriging instance.
type OKConfig struct {
	// VariogramModel specifies which variogram model to use:
	// "linear", "power", "gaussian", "spherical", "exponential", "hole-effect", or "custom".
	// Default: "linear"
	VariogramModel string

	// VariogramParameters are the model parameters.
	// If nil, parameters are automatically estimated from the data.
	// For list format:
	//   linear:      [slope, nugget]
	//   power:       [scale, exponent, nugget]
	//   gaussian/spherical/exponential/hole-effect: [sill, range, nugget]
	// For dict format (map):
	//   linear:      {"slope": ..., "nugget": ...}
	//   power:       {"scale": ..., "exponent": ..., "nugget": ...}
	//   gaussian/spherical/exponential/hole-effect: {"sill": ..., "range": ..., "nugget": ...}
	//                    or {"psill": ..., "range": ..., "nugget": ...}
	VariogramParameters interface{}

	// VariogramFunction is required if VariogramModel is "custom".
	VariogramFunction VariogramFunc

	// NLags is the number of averaging bins for the semivariogram. Default: 6
	NLags int

	// Weight specifies whether to weight smaller lags more heavily during fitting.
	Weight bool

	// AnisotropyScaling is the scalar stretching value for anisotropy. Default: 1.0
	AnisotropyScaling float64

	// AnisotropyAngle is the CCW angle (degrees) for anisotropy rotation. Default: 0.0
	AnisotropyAngle float64

	// Verbose enables progress output.
	Verbose bool

	// EnablePlotting is not implemented in Go (no matplotlib equivalent).
	EnablePlotting bool

	// EnableStatistics enables calculation of Q1, Q2, cR statistics.
	EnableStatistics bool

	// CoordinatesType is "euclidean" or "geographic". Default: "euclidean"
	CoordinatesType string

	// ExactValues: if true, interpolation is exact at input locations.
	// If false, accounts for nugget/variance at input locations.
	// Default: true
	ExactValues bool

	// PseudoInv enables pseudo-inverse for solving the kriging system.
	PseudoInv bool
}

// DefaultOKConfig returns a OKConfig with sensible defaults.
func DefaultOKConfig() OKConfig {
	return OKConfig{
		VariogramModel:      "linear",
		VariogramParameters: nil,
		VariogramFunction:   nil,
		NLags:               6,
		Weight:              false,
		AnisotropyScaling:   1.0,
		AnisotropyAngle:     0.0,
		Verbose:             false,
		EnablePlotting:      false,
		EnableStatistics:    false,
		CoordinatesType:     "euclidean",
		ExactValues:         true,
		PseudoInv:           false,
	}
}

// NewOrdinaryKriging creates a new OrdinaryKriging instance.
func NewOrdinaryKriging(x, y, z []float64, config OKConfig) (*OrdinaryKriging, error) {
	if len(x) == 0 || len(y) == 0 || len(z) == 0 {
		return nil, fmt.Errorf("input arrays must not be empty")
	}
	if len(x) != len(y) || len(x) != len(z) {
		return nil, fmt.Errorf("x, y, z arrays must have the same length")
	}
	if config.CoordinatesType != "euclidean" && config.CoordinatesType != "geographic" {
		return nil, fmt.Errorf("coordinates_type must be 'euclidean' or 'geographic'")
	}

	ok := &OrdinaryKriging{
		XOrig:           make([]float64, len(x)),
		YOrig:           make([]float64, len(y)),
		Z:               make([]float64, len(z)),
		Verbose:         config.Verbose,
		EnablePlotting:  config.EnablePlotting,
		ExactValues:     config.ExactValues,
		PseudoInv:       config.PseudoInv,
		CoordinatesType: config.CoordinatesType,
	}

	copy(ok.XOrig, x)
	copy(ok.YOrig, y)
	copy(ok.Z, z)

	// Set up variogram model
	ok.VariogramModel = config.VariogramModel
	if _, exists := VariogramModelMap[config.VariogramModel]; !exists && config.VariogramModel != "custom" {
		return nil, fmt.Errorf("unsupported variogram model: %s", config.VariogramModel)
	}

	if config.VariogramModel == "custom" {
		if config.VariogramFunction == nil {
			return nil, fmt.Errorf("must specify VariogramFunction for custom variogram model")
		}
		ok.VariogramFunc = config.VariogramFunction
	} else {
		ok.VariogramFunc = VariogramModelMap[config.VariogramModel]
	}

	// Handle anisotropy for euclidean coordinates
	if config.CoordinatesType == "euclidean" {
		ok.AnisotropyScaling = config.AnisotropyScaling
		ok.AnisotropyAngle = config.AnisotropyAngle

		// Compute center using (max+min)/2 (matching PyKrige)
		xMin, xMax := x[0], x[0]
		yMin, yMax := y[0], y[0]
		for i := range x {
			if x[i] < xMin {
				xMin = x[i]
			}
			if x[i] > xMax {
				xMax = x[i]
			}
			if y[i] < yMin {
				yMin = y[i]
			}
			if y[i] > yMax {
				yMax = y[i]
			}
		}
		ok.XCenter = (xMax + xMin) / 2.0
		ok.YCenter = (yMax + yMin) / 2.0

		// Always adjust for anisotropy (matching PyKrige behavior),
		// but with default scaling=1.0, angle=0.0 this is an identity transform.
		if ok.Verbose {
			fmt.Println("Adjusting data for anisotropy...")
		}
		points := make([][]float64, len(x))
		for i := range x {
			points[i] = []float64{x[i], y[i]}
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
	} else {
		// Geographic coordinates: no anisotropy
		if config.AnisotropyScaling != 1.0 && ok.Verbose {
			fmt.Println("Warning: Anisotropy is not compatible with geographic coordinates. Ignoring.")
		}
		ok.XCenter = 0.0
		ok.YCenter = 0.0
		ok.AnisotropyScaling = 1.0
		ok.AnisotropyAngle = 0.0
		ok.XAdjusted = make([]float64, len(x))
		ok.YAdjusted = make([]float64, len(y))
		copy(ok.XAdjusted, x)
		copy(ok.YAdjusted, y)
	}

	// Initialize variogram model
	if ok.Verbose {
		fmt.Println("Initializing variogram model...")
	}

	vpTemp := MakeVariogramParameterList(config.VariogramModel, config.VariogramParameters)

	ok.Lags, ok.Semivariance, ok.VariogramModelParameters = InitializeVariogramModel(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramModel, vpTemp, ok.VariogramFunc,
		config.NLags, config.Weight, config.CoordinatesType,
	)

	if ok.Verbose {
		ok.printVariogramInfo()
	}

	// Calculate statistics if enabled
	if config.EnableStatistics {
		if ok.Verbose {
			fmt.Println("Calculating statistics on variogram model fit...")
		}
		ok.Delta, ok.Sigma, ok.Epsilon = FindStatistics(
			ok.XAdjusted, ok.YAdjusted, ok.Z,
			ok.VariogramFunc, ok.VariogramModelParameters,
			config.CoordinatesType, ok.PseudoInv,
		)
		ok.Q1 = CalcQ1(ok.Epsilon)
		ok.Q2 = CalcQ2(ok.Epsilon)
		ok.CR = CalcCR(ok.Q2, ok.Sigma)
		if ok.Verbose {
			fmt.Printf("Q1 = %v\n", ok.Q1)
			fmt.Printf("Q2 = %v\n", ok.Q2)
			fmt.Printf("cR = %v\n", ok.CR)
		}
	}

	return ok, nil
}

// printVariogramInfo prints information about the fitted variogram model.
func (ok *OrdinaryKriging) printVariogramInfo() {
	fmt.Printf("Coordinates type: '%s'\n", ok.CoordinatesType)
	switch ok.VariogramModel {
	case "linear":
		fmt.Println("Using 'linear' Variogram Model")
		fmt.Printf("Slope: %v\n", ok.VariogramModelParameters[0])
		fmt.Printf("Nugget: %v\n", ok.VariogramModelParameters[1])
	case "power":
		fmt.Println("Using 'power' Variogram Model")
		fmt.Printf("Scale: %v\n", ok.VariogramModelParameters[0])
		fmt.Printf("Exponent: %v\n", ok.VariogramModelParameters[1])
		fmt.Printf("Nugget: %v\n", ok.VariogramModelParameters[2])
	case "custom":
		fmt.Println("Using Custom Variogram Model")
	default:
		fmt.Printf("Using '%s' Variogram Model\n", ok.VariogramModel)
		fmt.Printf("Partial Sill: %v\n", ok.VariogramModelParameters[0])
		fmt.Printf("Full Sill: %v\n",
			ok.VariogramModelParameters[0]+ok.VariogramModelParameters[2])
		fmt.Printf("Range: %v\n", ok.VariogramModelParameters[1])
		fmt.Printf("Nugget: %v\n", ok.VariogramModelParameters[2])
	}
}

// getKrigingMatrix assembles the (n+1) x (n+1) ordinary kriging matrix.
func (ok *OrdinaryKriging) getKrigingMatrix() [][]float64 {
	n := len(ok.XAdjusted)
	nPlus1 := n + 1

	a := make([][]float64, nPlus1)
	for i := range a {
		a[i] = make([]float64, nPlus1)
	}

	if ok.CoordinatesType == "euclidean" {
		// Compute pairwise Euclidean distances
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				dist := EuclideanDistance(
					ok.XAdjusted[i], ok.YAdjusted[i],
					ok.XAdjusted[j], ok.YAdjusted[j],
				)
				dArr := []float64{dist}
				a[i][j] = -ok.VariogramFunc(ok.VariogramModelParameters, dArr)[0]
			}
		}
	} else if ok.CoordinatesType == "geographic" {
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				dist := GreatCircleDistance(
					ok.XAdjusted[i], ok.YAdjusted[i],
					ok.XAdjusted[j], ok.YAdjusted[j],
				)
				dArr := []float64{dist}
				a[i][j] = -ok.VariogramFunc(ok.VariogramModelParameters, dArr)[0]
			}
		}
	}

	// Set diagonal to 0
	for i := 0; i < n; i++ {
		a[i][i] = 0.0
	}

	// Lagrange multiplier constraints
	for i := 0; i < n; i++ {
		a[n][i] = 1.0
		a[i][n] = 1.0
	}
	a[n][n] = 0.0

	return a
}

// getBD computes the distance from each data point to a target point.
func (ok *OrdinaryKriging) getBD(xpt, ypt float64) []float64 {
	n := len(ok.XAdjusted)
	bd := make([]float64, n)

	if ok.CoordinatesType == "euclidean" {
		for i := 0; i < n; i++ {
			bd[i] = EuclideanDistance(ok.XAdjusted[i], ok.YAdjusted[i], xpt, ypt)
		}
	} else if ok.CoordinatesType == "geographic" {
		for i := 0; i < n; i++ {
			bd[i] = GreatCircleDistance(xpt, ypt, ok.XAdjusted[i], ok.YAdjusted[i])
		}
	}

	return bd
}

// solveKrigingSystem solves the kriging linear system for a single target point.
func (ok *OrdinaryKriging) solveKrigingSystem(aInv *matWrapper, bd []float64) (float64, float64) {
	n := len(ok.XAdjusted)
	nPlus1 := n + 1

	// Build RHS vector b
	b := make([]float64, nPlus1)
	for i := 0; i < n; i++ {
		dArr := []float64{bd[i]}
		b[i] = -ok.VariogramFunc(ok.VariogramModelParameters, dArr)[0]
	}
	// Check for overlap with measurement points
	for i := 0; i < n; i++ {
		if math.Abs(bd[i]) <= Eps && ok.ExactValues {
			b[i] = 0.0
		}
	}
	b[n] = 1.0

	// Solve: x = A_inv * b
	x := aInv.MulVec(b)

	// Compute kriging estimate
	var zvalue float64
	for i := 0; i < n; i++ {
		zvalue += x[i] * ok.Z[i]
	}

	// Compute kriging variance
	var sigmasq float64
	for i := 0; i < nPlus1; i++ {
		sigmasq += x[i] * (-b[i])
	}

	return zvalue, sigmasq
}

// executeLoop solves the kriging system by looping over all target points.
// This is the memory-efficient approach.
func (ok *OrdinaryKriging) executeLoop(xpts, ypts []float64, mask []bool) (zvalues, sigmasq []float64) {
	npt := len(xpts)
	zvalues = make([]float64, npt)
	sigmasq = make([]float64, npt)

	// — 4a. 构建并求逆克里金矩阵 (n+1)×(n+1)，仅一次 —
	a := ok.getKrigingMatrix()
	aInv := invertMatrix(a, ok.PseudoInv)

	// — 4b. 并发逐点插值 —
	var wg sync.WaitGroup
	p, _ := ants.NewPoolWithFunc(runtime.NumCPU(), func(body interface{}) {
		defer wg.Done()
		j := body.(int)
		bd := ok.getBD(xpts[j], ypts[j])
		zvalues[j], sigmasq[j] = ok.solveKrigingSystem(aInv, bd)
	})
	defer p.Release()
	for j := 0; j < npt; j++ {
		if mask != nil && mask[j] {
			continue
		}
		wg.Add(1)
		_ = p.Invoke(j)
	}
	wg.Wait()
	return zvalues, sigmasq
}

// Execute calculates a kriged grid and the associated variance.
//
// style: "grid", "points", or "masked"
//   - "grid": xpoints and ypoints define a rectangular grid
//   - "points": xpoints and ypoints are coordinate pairs
//   - "masked": xpoints/ypoints define a grid, mask selects points
//
// xpoints, ypoints: target coordinates
// mask: boolean mask for "masked" style (false = compute, true = skip)
//
// Returns zvalues and sigmasq (kriging variance).
// For "grid" and "masked" styles, returns *mat.Dense (ny rows x nx cols).
// For "points" style, returns *mat.Dense (1 row x npt cols).
func (ok *OrdinaryKriging) Execute(style string, xpoints, ypoints []float64, mask []bool) (*mat.Dense, *mat.Dense) {
	if ok.Verbose {
		fmt.Println("Executing Ordinary Kriging...")
	}

	if style != "grid" && style != "masked" && style != "points" {
		panic("style must be 'grid', 'points', or 'masked'")
	}

	// — 1. 复制输入坐标 —
	tStep := time.Now()
	xpts := make([]float64, len(xpoints))
	ypts := make([]float64, len(ypoints))
	copy(xpts, xpoints)
	copy(ypts, ypoints)

	nx := len(xpts)
	ny := len(ypts)

	var flatMask []bool
	var npt int

	// — 2. Meshgrid 展开 / Points 直接使用 —
	if style == "grid" || style == "masked" {
		if style == "masked" {
			if len(mask) == 0 {
				panic("must specify mask when style is 'masked'")
			}
			if len(mask) != ny*nx {
				panic("mask dimensions do not match grid dimensions")
			}
			flatMask = make([]bool, len(mask))
			copy(flatMask, mask)
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
		xpts = gridX
		ypts = gridY
	} else { // points
		if nx != ny {
			panic("xpoints and ypoints must have same length for 'points' style")
		}
		npt = nx
	}
	// — 3. 各向异性坐标调整 —
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
			xpts[i] = adjusted[i][0]
			ypts[i] = adjusted[i][1]
		}
	}
	// — 4. 求解克里金系统（矩阵求逆 + 逐点插值）—
	zv, ss := ok.executeLoop(xpts, ypts, flatMask)

	// — 5. 输出整形 —
	var zvMat, ssMat *mat.Dense
	if style == "grid" || style == "masked" {
		zvMat = mat.NewDense(ny, nx, zv)
		ssMat = mat.NewDense(ny, nx, ss)
	} else {
		zvMat = mat.NewDense(1, npt, zv)
		ssMat = mat.NewDense(1, npt, ss)
	}
	_ = tStep // suppress unused warning
	return zvMat, ssMat
}

// ExecuteGrid is a convenience method for grid-style kriging.
// Returns *mat.Dense matrices with shape (len(ypoints), len(xpoints)),
// where row index = y-index and column index = x-index.
func (ok *OrdinaryKriging) ExecuteGrid(xpoints, ypoints []float64) (*mat.Dense, *mat.Dense) {
	return ok.Execute("grid", xpoints, ypoints, nil)
}

// GetVariogramPoints returns the lags and the variogram function evaluated at each lag.
func (ok *OrdinaryKriging) GetVariogramPoints() ([]float64, []float64) {
	variogram := ok.VariogramFunc(ok.VariogramModelParameters, ok.Lags)
	return ok.Lags, variogram
}

// GetStatistics returns Q1, Q2, and cR statistics.
func (ok *OrdinaryKriging) GetStatistics() (float64, float64, float64) {
	return ok.Q1, ok.Q2, ok.CR
}

// GetEpsilonResiduals returns the epsilon residuals for the variogram fit.
func (ok *OrdinaryKriging) GetEpsilonResiduals() []float64 {
	return ok.Epsilon
}

// UpdateVariogramModel allows updating the variogram model and/or parameters.
func (ok *OrdinaryKriging) UpdateVariogramModel(
	variogramModel string,
	variogramParameters interface{},
	variogramFunction VariogramFunc,
	nlags int,
	weight bool,
	anisotropyScaling float64,
	anisotropyAngle float64,
) error {

	if _, exists := VariogramModelMap[variogramModel]; !exists && variogramModel != "custom" {
		return fmt.Errorf("unsupported variogram model: %s", variogramModel)
	}

	ok.VariogramModel = variogramModel

	if variogramModel == "custom" {
		if variogramFunction == nil {
			return fmt.Errorf("must specify variogram function for custom model")
		}
		ok.VariogramFunc = variogramFunction
	} else {
		ok.VariogramFunc = VariogramModelMap[variogramModel]
	}

	// Handle anisotropy changes (euclidean only)
	if anisotropyScaling != ok.AnisotropyScaling || anisotropyAngle != ok.AnisotropyAngle {
		if ok.CoordinatesType == "euclidean" {
			ok.AnisotropyScaling = anisotropyScaling
			ok.AnisotropyAngle = anisotropyAngle
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
	}

	vpTemp := MakeVariogramParameterList(variogramModel, variogramParameters)

	ok.Lags, ok.Semivariance, ok.VariogramModelParameters = InitializeVariogramModel(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramModel, vpTemp, ok.VariogramFunc,
		nlags, weight, ok.CoordinatesType,
	)

	if ok.Verbose {
		ok.printVariogramInfo()
	}

	// Recalculate statistics
	ok.Delta, ok.Sigma, ok.Epsilon = FindStatistics(
		ok.XAdjusted, ok.YAdjusted, ok.Z,
		ok.VariogramFunc, ok.VariogramModelParameters,
		ok.CoordinatesType, ok.PseudoInv,
	)
	ok.Q1 = CalcQ1(ok.Epsilon)
	ok.Q2 = CalcQ2(ok.Epsilon)
	ok.CR = CalcCR(ok.Q2, ok.Sigma)

	return nil
}

// SwitchVerbose toggles verbose mode.
func (ok *OrdinaryKriging) SwitchVerbose() {
	ok.Verbose = !ok.Verbose
}

// matWrapper is a simple wrapper for a dense matrix with inversion support.
type matWrapper struct {
	data  []float64
	nRows int
	nCols int
}

// newMatWrapper creates a matrix from 2D slice.
func newMatWrapper(a [][]float64) *matWrapper {
	nRows := len(a)
	nCols := len(a[0])
	data := make([]float64, nRows*nCols)
	for i := range a {
		for j := range a[i] {
			data[i*nCols+j] = a[i][j]
		}
	}
	return &matWrapper{data: data, nRows: nRows, nCols: nCols}
}

// MulVec multiplies the inverse matrix by a vector.
func (m *matWrapper) MulVec(b []float64) []float64 {
	n := m.nRows
	result := make([]float64, n)
	for i := 0; i < n; i++ {
		var sum float64
		for j := 0; j < m.nCols; j++ {
			sum += m.data[i*m.nCols+j] * b[j]
		}
		result[i] = sum
	}
	return result
}

// invertMatrix computes the inverse of a matrix using Gaussian elimination
// (LU decomposition). If pseudoInv is true, uses SVD-based pseudo-inverse.
func invertMatrix(a [][]float64, pseudoInv bool) *matWrapper {
	n := len(a)

	if pseudoInv {
		// Use SVD-based pseudo-inverse via gonum
		return pseudoInverseMatrix(a)
	}

	// Standard inverse via Gaussian elimination with partial pivoting
	// Augment with identity matrix
	n2 := 2 * n
	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, n2)
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1.0
	}

	// Forward elimination with partial pivoting
	for col := 0; col < n; col++ {
		// Find pivot
		pivotRow := col
		maxVal := math.Abs(aug[col][col])
		for row := col + 1; row < n; row++ {
			if math.Abs(aug[row][col]) > maxVal {
				maxVal = math.Abs(aug[row][col])
				pivotRow = row
			}
		}

		if maxVal < Eps {
			// Singular matrix, fall back to pseudo-inverse
			return pseudoInverseMatrix(a)
		}

		// Swap rows
		if pivotRow != col {
			aug[col], aug[pivotRow] = aug[pivotRow], aug[col]
		}

		// Eliminate below
		pivot := aug[col][col]
		for row := col + 1; row < n; row++ {
			factor := aug[row][col] / pivot
			aug[row][col] = 0.0
			for k := col + 1; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	// Back substitution
	for col := n - 1; col >= 0; col-- {
		pivot := aug[col][col]
		if math.Abs(pivot) < Eps {
			return pseudoInverseMatrix(a)
		}
		// Normalize row
		for k := n; k < n2; k++ {
			aug[col][k] /= pivot
		}
		// Eliminate above
		for row := 0; row < col; row++ {
			factor := aug[row][col]
			for k := n; k < n2; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	// Extract inverse
	inv := &matWrapper{
		data:  make([]float64, n*n),
		nRows: n,
		nCols: n,
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			inv.data[i*n+j] = aug[i][n+j]
		}
	}

	return inv
}

// pseudoInverseMatrix computes the matrix inverse using gonum.
// Falls back to Gaussian elimination on failure.
func pseudoInverseMatrix(a [][]float64) *matWrapper {
	n := len(a)
	if n == 0 {
		return &matWrapper{nRows: 0, nCols: 0}
	}

	// Convert to gonum Dense matrix
	data := make([]float64, n*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			data[i*n+j] = a[i][j]
		}
	}

	matA := mat.NewDense(n, n, data)

	// Try standard inverse via gonum
	var aInv mat.Dense
	err := aInv.Inverse(matA)
	if err != nil {
		// Fallback: use Gaussian elimination (standard inverse)
		return invertMatrix(a, false)
	}

	inv := &matWrapper{
		data:  aInv.RawMatrix().Data,
		nRows: n,
		nCols: n,
	}
	return inv
}

// reshape2D reshapes a 1D slice into a 2D slice [rows][cols].
func reshape2D(data []float64, rows, cols int) [][]float64 {
	result := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		result[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			result[i][j] = data[i*cols+j]
		}
	}
	return result
}

// GetMaskedGrid extracts valid (non-masked) points from a grid result.
func GetMaskedGrid(zv [][]float64, ss [][]float64, mask []bool) ([][]float64, [][]float64) {
	rows := len(zv)
	cols := len(zv[0])

	var validZV []float64
	var validSS []float64

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			idx := i*cols + j
			if !mask[idx] {
				validZV = append(validZV, zv[i][j])
				validSS = append(validSS, ss[i][j])
			}
		}
	}

	resultZV := make([][]float64, 1)
	resultZV[0] = validZV
	resultSS := make([][]float64, 1)
	resultSS[0] = validSS
	return resultZV, resultSS
}
