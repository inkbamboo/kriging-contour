// Package kriging_contour 提供了基于克里金插值的等值线/等值面生成功能。
// 本文件包含了网格数据结构定义、克里金插值实现、网格操作等核心功能。
package kriging_contour

import (
	"encoding/json"
	"fmt"
	"math"
	"runtime"
	"sync"
	"time"

	"github.com/inkbamboo/kriging-contour/internal/kriging"
	"github.com/inkbamboo/kriging-contour/internal/utils"
	"github.com/panjf2000/ants/v2"
	"github.com/paulmach/orb"
	"gonum.org/v1/gonum/mat"
)

// Grid 表示克里金插值后的规则网格数据结构。
// 内部存储一个 *mat.Dense 矩阵（行对应 Y 坐标，列对应 X 坐标），
// 并保存对应的 X、Y 坐标轴。实现了 gonum/plot/plotter.GridXYZ 接口，
// 可直接用于等值线绘制和热力图渲染。
type Grid struct {
	Data    *mat.Dense // 插值后的网格数据矩阵（rows × cols）
	XCoords []float64  `json:"xCoords"` // X 坐标轴（列坐标）
	YCoords []float64  `json:"yCoords"` // Y 坐标轴（行坐标）
	Nx      int        `json:"nx"`      // X 方向的格点数
	Ny      int        `json:"ny"`      // Y 方向的格点数
}

// jsonGrid 是 Grid 的 JSON 序列化/反序列化中间结构。
// 将 mat.Dense 矩阵展开为一维 values 数组进行序列化。
type jsonGrid struct {
	Rows    int       `json:"rows"`
	Cols    int       `json:"cols"`
	Values  []float64 `json:"values"`
	XCoords []float64 `json:"xCoords"`
	YCoords []float64 `json:"yCoords"`
	Nx      int       `json:"nx"`
	Ny      int       `json:"ny"`
}

// MarshalJSON 将 Grid 序列化为 JSON。矩阵数据以行优先一维数组形式存储。
func (g *Grid) MarshalJSON() ([]byte, error) {
	rows, cols := g.Data.Dims()
	jg := jsonGrid{
		Rows:    rows,
		Cols:    cols,
		Values:  g.Data.RawMatrix().Data,
		XCoords: g.XCoords,
		YCoords: g.YCoords,
		Nx:      g.Nx,
		Ny:      g.Ny,
	}
	return json.Marshal(jg)
}

// UnmarshalJSON 从 JSON 反序列化 Grid。将一维数组还原为 mat.Dense 矩阵。
func (g *Grid) UnmarshalJSON(b []byte) error {
	var jg jsonGrid
	if err := json.Unmarshal(b, &jg); err != nil {
		return err
	}
	g.Data = mat.NewDense(jg.Rows, jg.Cols, jg.Values)
	g.XCoords = jg.XCoords
	g.YCoords = jg.YCoords
	g.Nx = jg.Nx
	g.Ny = jg.Ny
	return nil
}

// NewGrid 根据坐标轴和插值矩阵创建新的 Grid。
func NewGrid(xCoords, yCoords []float64, z *mat.Dense) *Grid {
	return &Grid{
		Data:    z,
		XCoords: xCoords,
		YCoords: yCoords,
		Nx:      len(xCoords),
		Ny:      len(yCoords),
	}
}

// generateGridCoords 根据边界范围和分辨率生成均匀的网格坐标轴。
// 在边界外增加 2% 的 padding，并按边界的长宽比例自动计算 nx 和 ny。
// resolution 指定较长边的格点数。
func generateGridCoords(minX, minY, maxX, maxY float64, resolution int) ([]float64, []float64) {
	padX := (maxX - minX) * 0.02
	padY := (maxY - minY) * 0.02
	minX -= padX
	maxX += padX
	minY -= padY
	maxY += padY

	dx, dy := maxX-minX, maxY-minY

	var nx, ny int
	if dx > dy {
		nx = resolution
		ny = int(math.Max(float64(resolution)*dy/dx, 10))
	} else {
		ny = resolution
		nx = int(math.Max(float64(resolution)*dx/dy, 10))
	}

	return utils.MakeLinespace(minX, maxX, nx), utils.MakeLinespace(minY, maxY, ny)
}

// Dims 返回网格的 X 和 Y 维度（列数, 行数）。实现 plotter.GridXYZ 接口。
func (g *Grid) Dims() (c, r int) { return g.Nx, g.Ny }

// Z 返回网格中 (c, r) 位置的值。实现 plotter.GridXYZ 接口。
func (g *Grid) Z(c, r int) float64 { return g.Data.At(r, c) }

// X 返回第 c 列的 X 坐标值。实现 plotter.GridXYZ 接口。
func (g *Grid) X(c int) float64 { return g.XCoords[c] }

// Y 返回第 r 行的 Y 坐标值。实现 plotter.GridXYZ 接口。
func (g *Grid) Y(r int) float64 { return g.YCoords[r] }

// newGridFromIDW 使用反距离加权（IDW）插值法生成网格。
// 当克里金插值失败时作为回退方案。使用 power=2 的权重函数，
// 对坐标进行归一化以平衡 X、Y 方向的影响。
func newGridFromIDW(px, py, pz, gridX, gridY []float64) *Grid {
	rows := len(gridY)
	cols := len(gridX)
	if rows == 0 || cols == 0 {
		return nil
	}

	minX, maxX := px[0], px[0]
	minY, maxY := py[0], py[0]
	for i := range px {
		if px[i] < minX {
			minX = px[i]
		}
		if px[i] > maxX {
			maxX = px[i]
		}
		if py[i] < minY {
			minY = py[i]
		}
		if py[i] > maxY {
			maxY = py[i]
		}
	}
	scaleX := maxX - minX
	scaleY := maxY - minY
	if scaleX < 1e-10 {
		scaleX = 1.0
	}
	if scaleY < 1e-10 {
		scaleY = 1.0
	}

	const power = 2.0
	zData := make([]float64, rows*cols)
	for r := 0; r < rows; r++ {
		for c := 0; c < cols; c++ {
			tx := gridX[c]
			ty := gridY[r]
			var weightSum, valueSum float64

			for i := range px {
				dx := (tx - px[i]) / scaleX
				dy := (ty - py[i]) / scaleY
				dist := math.Sqrt(dx*dx + dy*dy)

				if dist < 1e-12 {
					valueSum = pz[i]
					weightSum = 1.0
					break
				}

				w := 1.0 / math.Pow(dist, power)
				weightSum += w
				valueSum += w * pz[i]
			}

			if weightSum > 0 {
				zData[r*cols+c] = valueSum / weightSum
			} else {
				meanZ := 0.0
				for i := range pz {
					meanZ += pz[i]
				}
				zData[r*cols+c] = meanZ / float64(len(pz))
			}
		}
	}

	return &Grid{
		Data:    mat.NewDense(rows, cols, zData),
		XCoords: gridX,
		YCoords: gridY,
		Nx:      cols,
		Ny:      rows,
	}
}

// IdwFallBack 执行 IDW 回退插值并将结果归一化到原始数据的 Z 值范围。
// 当克里金插值失败或结果包含 NaN/Inf 时调用。
func IdwFallBack(x, y, z, gridX, gridY []float64, levels []float64) (grid *Grid, err error) {
	grid = newGridFromIDW(x, y, z, gridX, gridY)
	if grid == nil {
		return nil, fmt.Errorf("IDW fallback 插值也失败: 无法生成网格")
	}
	origMin, origMax := utils.MinMax(z...)
	gridMin, gridMax := utils.MinMax(grid.Data.RawMatrix().Data...)
	grid.Normalize(gridMin, gridMax, origMin, origMax, levels)
	return grid, nil
}

// Normalize 将网格值从 [gridMin, gridMax] 线性映射到目标范围。
// 优先使用 levels 的首尾值作为目标范围；若 levels 范围过小，
// 则使用 origMin/origMax。用于将克里金插值结果恢复到原始数据的 Z 值范围。
func (g *Grid) Normalize(gridMin, gridMax, origMin, origMax float64, levels []float64) {
	rows, cols := g.Data.Dims()
	if gridMax-gridMin <= 1e-10 {
		return
	}
	if len(levels) > 0 && levels[len(levels)-1]-levels[0] > 1e-10 {
		for i := 0; i < rows; i++ {
			for j := 0; j < cols; j++ {
				v := g.Data.At(i, j)
				v = levels[0] + (v-gridMin)/(gridMax-gridMin)*(levels[len(levels)-1]-levels[0])
				g.Data.Set(i, j, v)
			}
		}
	} else if origMax-origMin > 1e-10 {
		for i := 0; i < rows; i++ {
			for j := 0; j < cols; j++ {
				v := g.Data.At(i, j)
				v = origMin + (v-gridMin)/(gridMax-gridMin)*(origMax-origMin)
				g.Data.Set(i, j, v)
			}
		}
	}
}

// Sample 在网格的任意位置 (lon, lat) 进行双线性插值采样。
// 使用二分查找定位所在网格单元，然后进行双线性插值。
// 超出网格范围的点会被 clamp 到最近的边界单元格。
func (g *Grid) Sample(lon, lat float64) float64 {
	if g.Nx < 2 || g.Ny < 2 {
		return 0
	}

	j := sortSearch(g.Nx, func(j int) bool { return g.XCoords[j] >= lon })
	if j >= g.Nx {
		j = g.Nx - 1
	}
	if j > 0 && (j >= g.Nx || lon < g.XCoords[j]) {
		j--
	}
	if j >= g.Nx-1 {
		j = g.Nx - 2
	}
	if j < 0 {
		j = 0
	}

	i := sortSearch(g.Ny, func(i int) bool { return g.YCoords[i] >= lat })
	if i >= g.Ny {
		i = g.Ny - 1
	}
	if i > 0 && (i >= g.Ny || lat < g.YCoords[i]) {
		i--
	}
	if i >= g.Ny-1 {
		i = g.Ny - 2
	}
	if i < 0 {
		i = 0
	}

	x0, x1 := g.XCoords[j], g.XCoords[j+1]
	y0, y1 := g.YCoords[i], g.YCoords[i+1]
	if x1-x0 < 1e-12 || y1-y0 < 1e-12 {
		return g.Z(j, i)
	}

	tx := (lon - x0) / (x1 - x0)
	ty := (lat - y0) / (y1 - y0)
	tx = math.Max(0, math.Min(1, tx))
	ty = math.Max(0, math.Min(1, ty))

	v00 := g.Z(j, i)
	v10 := g.Z(j+1, i)
	v01 := g.Z(j, i+1)
	v11 := g.Z(j+1, i+1)

	return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
}

// sortSearch 二分查找，返回第一个使 f(i) 为 true 的 i。
// 功能等同于 sort.Search。
func sortSearch(n int, f func(int) bool) int {
	lo, hi := 0, n
	for lo < hi {
		mid := int(uint(lo+hi) >> 1)
		if !f(mid) {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

// GenerateGrid 执行克里金插值生成规则网格。
//
// 完整流程：
//  1. 根据边界多边形计算网格范围并生成坐标轴
//  2. 使用 DefaultOKConfig（spherical 模型，euclidean 坐标）创建 OrdinaryKriging
//  3. 执行 ExecuteGrid 获得插值矩阵
//  4. 若克里金插值失败或结果异常，回退到 IDW 插值
//  5. 对结果进行 Normalize，拉伸到原始 Z 值范围
//
// 注意：X/Y 坐标在传入克里金模块时发生了交换（经度→Y，纬度→X），
// 这是为了匹配 orb 库的 (lat, lon) 坐标约定。
func GenerateGrid(points []*Point, boundary orb.Polygon, resolution int, levels []float64) (grid *Grid, err error) {
	start := time.Now()
	fmt.Printf("[INFO] 步骤开始: 克里金网格插值\n")
	defer func() {
		if err != nil {
			fmt.Printf("[ERROR] 步骤失败: 克里金网格插值 (耗时: %v) - %v\n", time.Since(start).Round(time.Millisecond), err)
		} else {
			fmt.Printf("[INFO] 步骤完成: 克里金网格插值 (耗时: %v)\n", time.Since(start).Round(time.Millisecond))
		}
	}()
	minLon, maxLon := boundary.Bound().Min[1], boundary.Bound().Max[1]
	minLat, maxLat := boundary.Bound().Min[0], boundary.Bound().Max[0]
	gridX, gridY := generateGridCoords(minLon, minLat, maxLon, maxLat, resolution)
	x := make([]float64, len(points))
	y := make([]float64, len(points))
	z := make([]float64, len(points))
	for i, p := range points {
		x[i] = p.Y
		y[i] = p.X
		z[i] = p.Z
	}
	var ok *kriging.OrdinaryKriging
	config := kriging.DefaultOKConfig()
	config.CoordinatesType = "euclidean"
	config.VariogramModel = "spherical"
	config.Verbose = false
	ok, err = kriging.NewOrdinaryKriging(x, y, z, config)
	if err != nil || ok == nil {
		err = fmt.Errorf("克里金插值失败 (将使用 IDW fallback): %w", err)
		fmt.Printf("[WARN] 克里金插值失败，启用 IDW fallback 插值\n")
		return IdwFallBack(x, y, z, gridX, gridY, levels)
	}
	var zvGrid *mat.Dense
	zvGrid, _ = ok.ExecuteGrid(gridX, gridY)
	if zvGrid == nil {
		err = fmt.Errorf("克里金网格执行失败 (将使用 IDW fallback): %w", err)
		fmt.Printf("[WARN] 克里金网格执行失败，启用 IDW fallback 插值\n")
		return IdwFallBack(x, y, z, gridX, gridY, levels)
	}

	origMin, origMax := utils.MinMax(z...)
	gridMin, gridMax := utils.MinMax(zvGrid.RawMatrix().Data...)

	if math.IsNaN(gridMin) || math.IsNaN(gridMax) || math.IsInf(gridMin, 0) || math.IsInf(gridMax, 0) {
		fmt.Printf("[WARN] 克里金插值结果包含 NaN/Inf，使用 IDW fallback\n")
		return IdwFallBack(x, y, z, gridX, gridY, levels)
	}

	grid = NewGrid(gridX, gridY, zvGrid)
	grid.Normalize(gridMin, gridMax, origMin, origMax, levels)
	return grid, nil
}

// CombineGrid 将多个网格数据合并为单个网格。
//
// 通过逐点应用转换函数 convFn 和聚合函数 combineFn，生成新的网格数据。
// 支持并行处理以提高大网格的合并性能。
//
// 参数：
//   - gridList: 待合并的网格切片，所有网格必须具有相同的尺寸
//   - convFn: 转换函数，接收网格索引和该网格在当前位置的值，返回转换后的值
//   - combineFn: 聚合函数，接收每个网格转换后的值切片，返回合并结果
//
// 返回：
//   - *Grid: 合并后的网格，坐标轴与输入网格相同
//   - error: 错误信息，包括网格列表为空、网格为空或尺寸不一致等情况
//
// 实现细节：
//   1. 验证输入网格列表非空且所有网格具有相同尺寸
//   2. 创建结果网格矩阵，尺寸与输入网格相同
//   3. 使用 goroutine 池并行处理每一行数据
//   4. 对每个网格点，先应用 convFn 转换，再应用 combineFn 聚合
//   5. 使用预分配的切片避免重复内存分配
//   6. 使用互斥锁保护对结果矩阵的并发写入
func CombineGrid(gridList []*Grid, convFn func(gridIdx int, gridVal float64) float64, combineFn func(vals []float64) float64) (*Grid, error) {
	if len(gridList) == 0 {
		return nil, fmt.Errorf("gridList is empty")
	}
	refGrid := gridList[0]
	if refGrid == nil {
		return nil, fmt.Errorf("")
	}
	cols, rows := refGrid.Dims()

	// 检查所有网格尺寸是否一致
	for _, grid := range gridList {
		if grid == nil {
			return nil, fmt.Errorf("")
		}
		c, r := grid.Dims()
		if c != cols || r != rows {
			// 尺寸不一致，无法合并
			return nil, fmt.Errorf("")
		}
	}
	// 直接使用引用网格的坐标切片
	xCoords := refGrid.XCoords
	yCoords := refGrid.YCoords
	combinedData := mat.NewDense(rows, cols, nil)
	poolSize := utils.Min(runtime.NumCPU(), rows)
	if poolSize <= 0 {
		poolSize = 1
	}
	var wg sync.WaitGroup
	var lock sync.Mutex

	// 使用 goroutine 池并行处理每个 level 的等值线提取
	p, _ := ants.NewPoolWithFunc(poolSize, func(body interface{}) {
		defer wg.Done()
		// 为当前行预分配行数据
		rowData := make([]float64, cols)
		// 预分配 vals 切片，避免在循环中重复分配
		vals := make([]float64, len(gridList))
		rowIdx := body.(int)
		for c := 0; c < cols; c++ {
			// 收集所有网格在当前位置的值
			for idx, grid := range gridList {
				vals[idx] = convFn(idx, grid.Z(c, rowIdx))
			}
			rowData[c] = combineFn(vals)
		}
		// 使用互斥锁保护对 combinedData 的并发访问
		lock.Lock()
		combinedData.SetRow(rowIdx, rowData)
		lock.Unlock()
	})
	defer p.Release()
	for row := 0; row < rows; row++ {
		wg.Add(1)
		_ = p.Invoke(row)
	}
	wg.Wait()
	return &Grid{
		Data:    combinedData,
		XCoords: xCoords,
		YCoords: yCoords,
		Nx:      cols,
		Ny:      rows,
	}, nil
}
