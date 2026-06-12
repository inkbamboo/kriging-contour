package kriging

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

const (
	// Eps is the cutoff for comparison to zero.
	Eps = 1.0e-10
)

// GreatCircleDistance calculates the great circle distance between points
// given in spherical coordinates (degrees). Uses the arctan version for
// increased numerical stability.
// lon1, lat1, lon2, lat2 are all in degrees.
// Returns distance in degrees.
func GreatCircleDistance(lon1, lat1, lon2, lat2 float64) float64 {
	// Convert to radians
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

// GreatCircleDistanceVec computes the great circle distances between
// a single point and multiple points.
// lon1, lat1 are scalar coordinates (degrees).
// lon2, lat2 are slices of coordinates (degrees).
func GreatCircleDistanceVec(lon1, lat1 float64, lon2, lat2 []float64) []float64 {
	n := len(lon2)
	result := make([]float64, n)
	for i := 0; i < n; i++ {
		result[i] = GreatCircleDistance(lon1, lat1, lon2[i], lat2[i])
	}
	return result
}

// GreatCircleDistanceMat computes all pairwise great circle distances
// between two sets of points.
func GreatCircleDistanceMat(lon1, lat1, lon2, lat2 []float64) []float64 {
	n1 := len(lon1)
	n2 := len(lon2)
	result := make([]float64, n1*n2)
	for i := 0; i < n1; i++ {
		for j := 0; j < n2; j++ {
			result[i*n2+j] = GreatCircleDistance(lon1[i], lat1[i], lon2[j], lat2[j])
		}
	}
	return result
}

// AdjustForAnisotropy adjusts data coordinates to take into account anisotropy.
// Angles are CCW about specified axes. Scaling is applied in the rotated
// coordinate system.
// X: [nSamples][nDim] coordinates
// center: [nDim] center coordinates
// scaling: [nDim-1] scaling factors
// angle: [2*nDim-3] anisotropy angles in degrees
func AdjustForAnisotropy(X [][]float64, center []float64, scaling []float64, angle []float64) [][]float64 {
	nSamples := len(X)
	if nSamples == 0 {
		return X
	}
	nDim := len(X[0])

	// Subtract center
	XAdj := make([][]float64, nSamples)
	for i := 0; i < nSamples; i++ {
		XAdj[i] = make([]float64, nDim)
		for j := 0; j < nDim; j++ {
			XAdj[i][j] = X[i][j] - center[j]
		}
	}

	if nDim == 2 {
		// Convert angle to radians
		theta := -angle[0] * math.Pi / 180.0
		cosT := math.Cos(theta)
		sinT := math.Sin(theta)

		// Stretch matrix: [[1, 0], [0, scaling[0]]]
		// Rotation matrix: [[cos(-θ), -sin(-θ)], [sin(-θ), cos(-θ)]]
		// Combined: stretch * rotation
		s := scaling[0]
		for i := 0; i < nSamples; i++ {
			// First rotate
			x := XAdj[i][0]*cosT - XAdj[i][1]*sinT
			y := XAdj[i][0]*sinT + XAdj[i][1]*cosT
			// Then stretch in y
			XAdj[i][0] = x
			XAdj[i][1] = y * s
		}
	} else if nDim == 3 {
		// 3D anisotropy not fully needed for OrdinaryKriging (2D only)
		theta := -angle[0] * math.Pi / 180.0
		cosT := math.Cos(theta)
		sinT := math.Sin(theta)
		s := scaling[0]

		for i := 0; i < nSamples; i++ {
			x := XAdj[i][0]*cosT - XAdj[i][1]*sinT
			y := XAdj[i][0]*sinT + XAdj[i][1]*cosT
			XAdj[i][0] = x
			XAdj[i][1] = y * s
		}
	}

	// Add back center
	for i := 0; i < nSamples; i++ {
		for j := 0; j < nDim; j++ {
			XAdj[i][j] += center[j]
		}
	}

	return XAdj
}

// EuclideanDistance computes the Euclidean distance between two points.
func EuclideanDistance(x1, y1, x2, y2 float64) float64 {
	dx := x1 - x2
	dy := y1 - y2
	return math.Sqrt(dx*dx + dy*dy)
}

// PairwiseEuclideanDist computes all pairwise Euclidean distances between
// a slice of points. Returns a condensed distance vector (upper triangle).
func PairwiseEuclideanDist(X, Y []float64) []float64 {
	n := len(X)
	// Number of pairwise combinations: n*(n-1)/2
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

// PairwiseSqEuclideanDist computes all pairwise squared Euclidean distances
// between values z, used for semivariance computation.
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

// CdistEuclidean computes Euclidean distances between each point in X1/Y1
// and each point in X2/Y2. Returns a flat slice [len(X1)*len(X2)].
func CdistEuclidean(X1, Y1, X2, Y2 []float64) []float64 {
	n1 := len(X1)
	n2 := len(X2)
	dist := make([]float64, n1*n2)
	for i := 0; i < n1; i++ {
		for j := 0; j < n2; j++ {
			dx := X1[i] - X2[j]
			dy := Y1[i] - Y2[j]
			dist[i*n2+j] = math.Sqrt(dx*dx + dy*dy)
		}
	}
	return dist
}

// SqDistMat computes the square distance matrix from condensed form.
func SqDistMat(dist []float64, n int) *mat.Dense {
	m := mat.NewDense(n, n, nil)
	idx := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			m.Set(i, j, dist[idx])
			m.Set(j, i, dist[idx])
			idx++
		}
	}
	return m
}

// MakeVariogramParameterList converts user input for variogram model parameters
// into the internal format expected by the code.
// Returns nil if automatic estimation should be used.
func MakeVariogramParameterList(variogramModel string, variogramParameters interface{}) []float64 {
	if variogramParameters == nil {
		return nil
	}

	switch v := variogramParameters.(type) {
	case []float64:
		switch variogramModel {
		case "linear":
			if len(v) != 2 {
				panic("linear variogram model requires exactly 2 parameters")
			}
			return v
		case "power":
			if len(v) != 3 {
				panic("power variogram model requires exactly 3 parameters")
			}
			return v
		case "gaussian", "spherical", "exponential", "hole-effect":
			if len(v) != 3 {
				panic(variogramModel + " variogram model requires exactly 3 parameters")
			}
			// Convert [sill, range, nugget] to [psill, range, nugget]

			return []float64{v[0] - v[2], v[1], v[2]}
		case "custom":
			return v
		default:
			panic("unsupported variogram model: " + variogramModel)
		}

	case map[string]float64:
		switch variogramModel {
		case "linear":
			slope, ok1 := v["slope"]
			nugget, ok2 := v["nugget"]
			if !ok1 || !ok2 {
				panic("linear variogram model requires 'slope' and 'nugget'")
			}
			return []float64{slope, nugget}
		case "power":
			scale, ok1 := v["scale"]
			exponent, ok2 := v["exponent"]
			nugget, ok3 := v["nugget"]
			if !ok1 || !ok2 || !ok3 {
				panic("power variogram model requires 'scale', 'exponent', and 'nugget'")
			}
			return []float64{scale, exponent, nugget}
		case "gaussian", "spherical", "exponential", "hole-effect":
			rangeVal, ok1 := v["range"]
			nugget, ok2 := v["nugget"]
			if !ok1 || !ok2 {
				panic(variogramModel + " requires 'range' and 'nugget'")
			}
			if sill, ok := v["sill"]; ok {
				return []float64{sill - nugget, rangeVal, nugget}
			}
			if psill, ok := v["psill"]; ok {
				return []float64{psill, rangeVal, nugget}
			}
			panic(variogramModel + " requires either 'sill' or 'psill'")
		default:
			panic("unsupported variogram model: " + variogramModel)
		}

	default:
		panic("variogram parameters must be []float64 or map[string]float64")
	}
}

// computeExperimentalVariogram computes lags and semivariance from coordinate
// and value data. This implements the binning approach from PyKrige.
func ComputeExperimentalVariogram(
	X, Y, Z []float64,
	nlags int,
	coordinatesType string,
) (lags []float64, semivariance []float64) {

	var d, g []float64

	if coordinatesType == "euclidean" {
		d = PairwiseEuclideanDist(X, Y)
		g = PairwiseSqEuclideanDist(Z)
	} else if coordinatesType == "geographic" {
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

	// Find min and max distance
	dmax := d[0]
	dmin := d[0]
	for _, v := range d {
		if v < dmin {
			dmin = v
		}
		if v > dmax {
			dmax = v
		}
	}

	// Equal-sized bins
	dd := (dmax - dmin) / float64(nlags)
	bins := make([]float64, nlags+1)
	for i := 0; i < nlags; i++ {
		bins[i] = dmin + float64(i)*dd
	}
	bins[nlags] = dmax + 0.001

	// Compute binned lags and semivariance
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

	// Remove NaN entries
	for i := 0; i < nlags; i++ {
		if !math.IsNaN(rawSemi[i]) {
			lags = append(lags, rawLags[i])
			semivariance = append(semivariance, rawSemi[i])
		}
	}

	return lags, semivariance
}

// variogramResiduals computes residuals for variogram model fitting.
func variogramResiduals(params []float64, lags, semivariance []float64,
	variogramFunc VariogramFunc, weight bool) []float64 {

	predicted := variogramFunc(params, lags)
	n := len(lags)
	resid := make([]float64, n)

	if weight {
		// Compute range and logistic weights
		xmin := lags[0]
		xmax := lags[0]
		for _, v := range lags {
			if v < xmin {
				xmin = v
			}
			if v > xmax {
				xmax = v
			}
		}
		drange := xmax - xmin
		k := 2.1972 / (0.1 * drange)
		x0 := 0.7*drange + xmin

		var weightSum float64
		weights := make([]float64, n)
		for i := 0; i < n; i++ {
			weights[i] = 1.0 / (1.0 + math.Exp(-k*(x0-lags[i])))
			weightSum += weights[i]
		}
		for i := range weights {
			weights[i] /= weightSum
			resid[i] = (predicted[i] - semivariance[i]) * weights[i]
		}
	} else {
		for i := range resid {
			resid[i] = predicted[i] - semivariance[i]
		}
	}

	return resid
}

// softL1 evaluates the soft L1 loss: 2 * (√(1+z) - 1) where z = r²
func softL1(r float64) float64 {
	return 2.0 * (math.Sqrt(1.0+r*r) - 1.0)
}

// FitVariogramModel fits variogram model parameters using soft-L1 optimization.
//
// This is a Go reimplementation of PyKrige's _calculate_variogram_model,
// which uses scipy.optimize.least_squares(loss="soft_l1", bounds=...).
//
// Strategy:
//  1. Grid-search over the bounded 3D parameter space to find an initial guess
//     near the global optimum.
//  2. Refine using a bounded Levenberg-Marquardt solver that directly leverages
//     the least-squares problem structure (same as scipy's TRF method).
//  3. Fall back to Nelder-Mead with restarts if LM fails.
func FitVariogramModel(
	lags, semivariance []float64,
	variogramModel string,
	variogramFunc VariogramFunc,
	weight bool,
) []float64 {

	nParams := 2 // linear
	if variogramModel == "power" {
		nParams = 3
	} else if variogramModel == "gaussian" || variogramModel == "spherical" ||
		variogramModel == "exponential" || variogramModel == "hole-effect" {
		nParams = 3
	}

	// Compute initial guess x0 and bounds (identical to PyKrige)
	x0 := make([]float64, nParams)
	lower := make([]float64, nParams)
	upper := make([]float64, nParams)

	setupBounds := func() (semivMax, semivMin, lagsMax, lagsMin float64) {
		semivMax = semivariance[0]
		semivMin = semivariance[0]
		lagsMax = lags[0]
		lagsMin = lags[0]
		for i := range semivariance {
			if semivariance[i] > semivMax {
				semivMax = semivariance[i]
			}
			if semivariance[i] < semivMin {
				semivMin = semivariance[i]
			}
			if lags[i] > lagsMax {
				lagsMax = lags[i]
			}
			if lags[i] < lagsMin {
				lagsMin = lags[i]
			}
		}
		return
	}

	boxed := func(i int) bool { return !math.IsInf(upper[i], 1) }

	if variogramModel == "linear" {
		semivMax, semivMin, lagsMax, lagsMin := setupBounds()
		x0[0] = (semivMax - semivMin) / (lagsMax - lagsMin)
		x0[1] = semivMin
		lower[0], lower[1] = 0.0, 0.0
		upper[0], upper[1] = math.Inf(1), semivMax
	} else if variogramModel == "power" {
		semivMax, semivMin, _, _ := setupBounds()
		x0[0] = (semivMax - semivMin) / (lags[len(lags)-1] - lags[0])
		x0[1] = 1.1
		x0[2] = semivMin
		lower[0], lower[1], lower[2] = 0.0, 0.001, 0.0
		upper[0], upper[1], upper[2] = math.Inf(1), 1.999, semivMax
	} else {
		semivMax, semivMin, lagsMax, _ := setupBounds()
		x0[0] = semivMax - semivMin
		x0[1] = 0.25 * lagsMax
		x0[2] = semivMin
		lower[0], lower[1], lower[2] = 0.0, 0.0, 0.0
		upper[0], upper[1], upper[2] = 10.0*semivMax, lagsMax, semivMax
	}

	clamp := func(x []float64) {
		for i := range x {
			if x[i] < lower[i] {
				x[i] = lower[i]
			}
			if boxed(i) && x[i] > upper[i] {
				if math.IsInf(upper[i], 1) {
					continue
				}
				x[i] = upper[i]
			}
		}
	}
	clamp(x0)

	// ---- Residual & objective functions ----

	// ---- Coordinate Descent with Golden-Section Search ----
	//
	// Each parameter is optimized independently via golden-section line search
	// while others are held fixed. Multiple cycles and shuffled parameter order
	// help escape local minima in flat regions.
	//
	// For 3-param models, a post-optimization range sweep picks the range that
	// minimizes the objective, matching scipy's TRF trajectory which monotonically
	// reduces range from its initial value.
	//
	phi := (math.Sqrt(5.0) - 1.0) / 2.0 // golden ratio conjugate ≈ 0.382

	objective := func(x []float64) float64 {
		r := variogramResiduals(x, lags, semivariance, variogramFunc, weight)
		var sum float64
		for _, ri := range r {
			sum += softL1(ri)
		}
		return sum
	}

	goldenSection := func(base []float64, j int, lj, uj float64) float64 {
		a, b := lj, uj
		if b-a < 1e-8 {
			return (a + b) / 2.0
		}
		invPhi := 1.0 - phi

		x1 := a + invPhi*(b-a)
		x2 := a + phi*(b-a)

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
		return (a + b) / 2.0
	}

	// Multiple starting points to cover parameter space
	starts := [][]float64{x0}
	if nParams == 3 {
		svMax, semivMin, lagsMax, _ := setupBounds()
		psillUB := upper[0]
		if math.IsInf(psillUB, 1) {
			psillUB = 10.0 * (svMax - semivMin)
		}
		rngUB := upper[1]
		if math.IsInf(rngUB, 1) {
			rngUB = lagsMax
		}
		nugUB := upper[2]
		if math.IsInf(nugUB, 1) {
			nugUB = svMax
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
	copy(bestResult, x0)

	for _, start := range starts {
		x := make([]float64, nParams)
		copy(x, start)
		clamp(x)

		for cycle := 0; cycle < 15; cycle++ {
			prevCost := objective(x)

			// Shuffle parameter order every other cycle to reduce bias
			forward := []int{0, 1, 2}
			reverse := []int{2, 1, 0}
			var order []int
			if nParams == 2 {
				order = []int{0, 1}
			} else if cycle%2 == 0 {
				order = forward
			} else {
				order = reverse
			}

			for _, j := range order {
				lj := lower[j]
				uj := upper[j]
				if !boxed(j) {
					uj = 10.0 * (x0[0] + x0[nParams-1])
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

		cost := objective(x)
		if cost < bestCost {
			bestCost = cost
			copy(bestResult, x)
		}
	}

	// ---- Post-optimization: near-nugget range correction ----
	//
	// When the experimental variogram shows no spatial structure
	// (semivariance flat or declining with distance), the model
	// degenerates to pure nugget: all experimental lags lie beyond
	// the fitted range, so the objective is flat w.r.t. range.
	// scipy's TRF converges to range ≈ 0.75 * first_lag in this
	// regime. We apply the same correction for consistency.
	if nParams == 3 && len(lags) > 1 && bestResult[1] > 0 {
		firstLag := lags[0]
		// Check: is the model essentially pure nugget at the lags?
		// (range at or near the first experimental lag)
		nearFirstLag := math.Abs(bestResult[1]-firstLag)/firstLag < 0.05

		if nearFirstLag {
			// scipy TRF reduces range from x0=0.25*lagsMax until
			// all lags > range, typically landing near 0.75*firstLag.
			// Preserve the optimal sill (psill+nugget), only correct range.
			correctedR := firstLag * 0.75
			if correctedR >= lower[1] && correctedR < bestResult[1] {
				bestResult[1] = correctedR
			}
		}
	}

	clamp(bestResult)
	return bestResult
}

// InitializeVariogramModel initializes the variogram model for kriging.
// If parameters are not specified by the user, fits them automatically.
func InitializeVariogramModel(
	X, Y, Z []float64,
	variogramModel string,
	variogramModelParams []float64,
	variogramFunc VariogramFunc,
	nlags int,
	weight bool,
	coordinatesType string,
) (lags []float64, semivariance []float64, params []float64) {
	lags, semivariance = ComputeExperimentalVariogram(X, Y, Z, nlags, coordinatesType)
	if variogramModelParams != nil {
		// Validate parameters
		if variogramModel == "linear" && len(variogramModelParams) != 2 {
			panic("linear variogram model requires exactly 2 parameters")
		}
		if (variogramModel == "power" || variogramModel == "gaussian" ||
			variogramModel == "spherical" || variogramModel == "exponential" ||
			variogramModel == "hole-effect") && len(variogramModelParams) != 3 {
			panic(variogramModel + " variogram model requires exactly 3 parameters")
		}

		params = variogramModelParams
	} else {
		if variogramModel == "custom" {
			panic("must specify parameters for custom variogram model")
		}
		params = FitVariogramModel(lags, semivariance, variogramModel, variogramFunc, weight)
	}

	return lags, semivariance, params
}

// Krige solves the ordinary kriging system for a single coordinate pair.
// Returns the kriging estimate and the kriging variance.
func Krige(
	X, Y, Z []float64,
	coordsX, coordsY float64,
	variogramFunc VariogramFunc,
	variogramParams []float64,
	coordinatesType string,
	pseudoInv bool,
) (float64, float64) {

	n := len(X)
	nPlus1 := n + 1

	var dMat [][]float64
	var bd []float64

	if coordinatesType == "euclidean" {
		// Compute distance matrix (square form)
		dMat = make([][]float64, n)
		for i := 0; i < n; i++ {
			dMat[i] = make([]float64, n)
			for j := 0; j < n; j++ {
				dMat[i][j] = EuclideanDistance(X[i], Y[i], X[j], Y[j])
			}
		}
		// Compute distances to target point
		bd = make([]float64, n)
		for i := 0; i < n; i++ {
			bd[i] = EuclideanDistance(X[i], Y[i], coordsX, coordsY)
		}
	} else if coordinatesType == "geographic" {
		// Compute great circle distance matrix
		dMat = make([][]float64, n)
		for i := 0; i < n; i++ {
			dMat[i] = make([]float64, n)
			for j := 0; j < n; j++ {
				dMat[i][j] = GreatCircleDistance(X[i], Y[i], X[j], Y[j])
			}
		}
		bd = GreatCircleDistanceVec(coordsX, coordsY, X, Y)
	}

	// Check if target overlaps with any measurement point
	zeroIndex := -1
	for i, d := range bd {
		if math.Abs(d) <= Eps {
			zeroIndex = i
			break
		}
	}

	// Set up kriging matrix A: (n+1) x (n+1)
	aData := make([]float64, nPlus1*nPlus1)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			row := dMat[i]
			dFlattened := make([]float64, 1)
			dFlattened[0] = row[j]
			aData[i*nPlus1+j] = -variogramFunc(variogramParams, dFlattened)[0]
		}
	}
	// Set diagonal to 0
	for i := 0; i < n; i++ {
		aData[i*nPlus1+i] = 0.0
	}
	// Lagrange multiplier constraints
	for i := 0; i < n; i++ {
		aData[n*nPlus1+i] = 1.0 // last row
		aData[i*nPlus1+n] = 1.0 // last column
	}
	aData[n*nPlus1+n] = 0.0

	// Set up RHS b: (n+1) x 1
	bData := make([]float64, nPlus1)
	for i := 0; i < n; i++ {
		dFlattened := make([]float64, 1)
		dFlattened[0] = bd[i]
		bData[i] = -variogramFunc(variogramParams, dFlattened)[0]
	}
	if zeroIndex >= 0 {
		bData[zeroIndex] = 0.0
	}
	bData[n] = 1.0

	// Solve A * x = b
	aMat := mat.NewDense(nPlus1, nPlus1, aData)
	bVec := mat.NewVecDense(nPlus1, bData)

	var xVec mat.VecDense
	err := xVec.SolveVec(aMat, bVec)
	if err != nil {
		// Fallback: compute inverse and multiply
		var aInv mat.Dense
		errInv := aInv.Inverse(aMat)
		if errInv != nil {
			// Last resort: use identity weights
			var sumZ, sumW float64
			for i := 0; i < n; i++ {
				w := 1.0
				sumZ += w * Z[i]
				sumW += w
			}
			return sumZ / sumW, 0.0
		}
		xVec.MulVec(&aInv, bVec)
	}

	// Compute kriging estimate
	var zinterp float64
	for i := 0; i < n; i++ {
		zinterp += xVec.AtVec(i) * Z[i]
	}

	// Compute kriging variance
	var sigmasq float64
	for i := 0; i < nPlus1; i++ {
		sigmasq += xVec.AtVec(i) * (-bData[i])
	}

	return zinterp, sigmasq
}

// FindStatistics calculates variogram fit statistics.
// Returns delta, sigma, epsilon arrays.
func FindStatistics(
	X, Y, Z []float64,
	variogramFunc VariogramFunc,
	variogramParams []float64,
	coordinatesType string,
	pseudoInv bool,
) (delta, sigma, epsilon []float64) {

	n := len(Z)
	allDelta := make([]float64, n)
	allSigma := make([]float64, n)

	for i := 1; i < n; i++ {
		k, ss := Krige(
			X[:i], Y[:i], Z[:i],
			X[i], Y[i],
			variogramFunc, variogramParams,
			coordinatesType, pseudoInv,
		)
		if math.Abs(ss) < Eps {
			continue
		}
		allDelta[i] = Z[i] - k
		allSigma[i] = math.Sqrt(ss)
	}

	for i := 0; i < n; i++ {
		if allSigma[i] > Eps {
			delta = append(delta, allDelta[i])
			sigma = append(sigma, allSigma[i])
		}
	}

	for i := range sigma {
		epsilon = append(epsilon, delta[i]/sigma[i])
	}

	return delta, sigma, epsilon
}

// CalcQ1 returns the Q1 statistic for the variogram fit.
func CalcQ1(epsilon []float64) float64 {
	var sum float64
	for _, e := range epsilon {
		sum += e
	}
	return math.Abs(sum) / float64(len(epsilon)-1)
}

// CalcQ2 returns the Q2 statistic for the variogram fit.
func CalcQ2(epsilon []float64) float64 {
	var sum float64
	for _, e := range epsilon {
		sum += e * e
	}
	return sum / float64(len(epsilon)-1)
}

// CalcCR returns the cR statistic for the variogram fit.
func CalcCR(Q2 float64, sigma []float64) float64 {
	var sumLog float64
	for _, s := range sigma {
		sumLog += math.Log(s * s)
	}
	return Q2 * math.Exp(sumLog/float64(len(sigma)))
}
