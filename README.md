# Kriging Contour - 克里金等值线/等值面生成库

一个基于克里金插值（Kriging）的等值线（Contour Lines）和等值面（Contour Faces）生成库，支持从离散采样点数据生成高质量的 GeoJSON 格式等值线和等值面。

## 功能特性

### 核心功能
- **等值线生成** (`GenerateLines`)：从离散点数据生成等值线，输出带 level、闭合状态、高值侧属性的 GeoJSON Feature
- **等值面生成** (`GenerateFaces`)：从等值线递归切割边界多边形，生成填充多边形等值面
- **克里金插值** (`GenerateGrid`)：使用 Ordinary Kriging 进行空间插值，支持多种变差函数模型
- **IDW 回退** (`IdwFallBack`)：当克里金插值失败时自动回退到反距离加权插值
- **地形特征提取** (`ExtractPeaks`)：基于 8 邻域比较算法检测山顶点（局部极大值）
- **网格合并** (`CombineGrid`)：通过自定义转换和聚合函数合并多个网格

### 数据处理
- **数据清洗** (`ValidatePoints`)：自动过滤 nil、NaN、Inf 和重复数据点
- **参数验证** (`ValidateContourOption`)：校验分辨率范围、LevelList/ContourInterval 配置
- **边界校验** (`ValidateBoundary`)：检测多边形自交、面积过小、顶点数不足
- **统计信息** (`ComputeDataStats`)：计算数据范围、均值、标准差、共线性等统计量

### 输出格式
- **GeoJSON 输出**：标准 GeoJSON Feature 格式，包含 level、is_closed、high_value_side 等属性
- **图像渲染** (`renderContourImage`)：可选的 1600×1200 像素等值线热力图 PNG 输出
- **网格序列化**：支持 Grid 的 JSON 序列化/反序列化

### 性能优化
- **并行处理**：使用 goroutine 池（ants）并行处理每个 level 的等值线提取和等值面 level 计算
- **空间哈希**：使用 cell grid 加速线段端点合并的邻近搜索
- **并查集**：Union-Find 算法合并连通分量
- **矩阵优化**：支持伪逆和高斯消元两种矩阵求逆策略

## 项目结构

```
kriging-contour/
├── type.go              # 数据结构定义（Point、ContourOption）
├── contour.go           # 等值线生成主逻辑（GenerateLines、线段合并、裁剪、定向）
├── faces.go             # 等值面生成算法（GenerateFaces、多边形切割、去重、level 计算）
├── grid.go              # 网格数据结构（Grid）和克里金插值（GenerateGrid、CombineGrid）
├── validate.go          # 数据校验和参数验证（ValidatePoints、ValidateBoundary 等）
├── peaks.go             # 地形特征提取（ExtractPeaks，8 邻域局部极大值检测）
├── image.go             # 等值线图像渲染（renderContourImage，热力图 PNG 输出）
├── canvas.go            # 自定义 vg.Canvas 实现（捕获 gonum/plot 的等值线路径）
├── go.mod               # Go 模块依赖
└── internal/            # 内部工具库
    ├── kriging/         # 克里金插值核心实现
    │   ├── core.go      # 距离计算、实验变差函数、变差函数拟合、单点克里金、交叉验证
    │   ├── ok.go        # Ordinary Kriging 实现（矩阵构建、执行插值、矩阵求逆）
    │   └── variogram.go # 变差函数模型定义（线性、幂、高斯、指数、球状、孔穴效应）
    └── utils/           # 通用工具函数
        └── util.go      # CloseRing、EqPoint、MakeLinespace、MinMax、MatToFlat
```

## 快速开始

### 安装

```bash
go get github.com/inkbamboo/kriging-contour
```

### 基本用法

```go
package main

import (
    "encoding/json"
    "fmt"
    "github.com/inkbamboo/kriging-contour"
    "github.com/paulmach/orb"
)

func main() {
    // 1. 准备输入数据
    points := []*kriging_contour.Point{
        {X: 120.0, Y: 30.0, Z: 100.0},
        {X: 120.1, Y: 30.1, Z: 150.0},
        // ... 更多数据点
    }
    
    // 2. 定义边界多边形
    boundary := orb.Polygon{{
        {120.0, 30.0}, {120.2, 30.0}, {120.2, 30.2}, {120.0, 30.2}, {120.0, 30.0},
    }}
    
    // 3. 配置参数
    opt := &kriging_contour.ContourOption{
        ContourType:      "temperature",
        ContourStart:     100.0,
        ContourEnd:       200.0,
        ContourInterval:  10.0,
        Resolution:       200,
        ImagePath:        "output/contour.png", // 可选
        Extra:            map[string]interface{}{"source": "weather_station"},
    }
    
    // 4. 生成等值线
    features := kriging_contour.GenerateLines(points, boundary, opt)
    
    // 5. 输出 GeoJSON
    geojson, _ := json.Marshal(features)
    fmt.Println(string(geojson))
}
```

### 等值面生成

```go
// 基于等值线生成等值面
features := kriging_contour.GenerateLines(points, boundary, opt)
faces := kriging_contour.GenerateFaces(features, boundary, opt)
```

### 基于已有网格生成等值线

```go
// 复用已有网格，跳过插值步骤
grid, _ := kriging_contour.GenerateGrid(points, boundary, 200, opt.LevelList)
features := kriging_contour.GenerateLinesByGrid(grid, boundary, opt)
```

### 提取山顶点

```go
grid, _ := kriging_contour.GenerateGrid(points, boundary, 200, opt.LevelList)
peaks := kriging_contour.ExtractPeaks(grid, boundary)
// peaks 按高度升序排列
```

### 合并多个网格

```go
// 取多个网格在每个位置的最大值
merged, err := kriging_contour.CombineGrid(
    gridList,
    func(gridIdx int, gridVal float64) float64 { return gridVal },
    func(vals []float64) float64 {
        maxVal := vals[0]
        for _, v := range vals[1:] {
            if v > maxVal { maxVal = v }
        }
        return maxVal
    },
)
```

## API 文档

### 主要类型

#### Point
空间采样数据点，包含 X（经度）、Y（纬度）、Z（观测值）坐标。

```go
type Point struct {
    X float64 `json:"x"` // 经度或 X 坐标
    Y float64 `json:"y"` // 纬度或 Y 坐标
    Z float64 `json:"z"` // 观测值
}
```

#### ContourOption
等值线/等值面生成参数配置。

```go
type ContourOption struct {
    ContourType      string                 `json:"contour_type"`          // 等值线类型名称
    ContourStart     float64                `json:"contour_start"`         // 起始数值
    ContourEnd       float64                `json:"contour_end,omitempty"` // 结束数值
    ContourInterval  float64                `json:"contour_interval"`      // 数值间隔
    LevelList        []float64              `json:"level_list,omitempty"`  // 显式级别列表（优先）
    Resolution       int                    `json:"resolution"`            // 网格分辨率（长边格点数）
    ImagePath        string                 `json:"image_path,omitempty"`  // 图片输出路径
    Extra            map[string]interface{} `json:"extra"`                 // 附加属性
}
```

#### Grid
克里金插值后的规则网格数据结构，实现 `plotter.GridXYZ` 接口。

```go
type Grid struct {
    Data    *mat.Dense // 插值后的网格数据矩阵（rows × cols）
    XCoords []float64  `json:"xCoords"` // X 坐标轴（列坐标）
    YCoords []float64  `json:"yCoords"` // Y 坐标轴（行坐标）
    Nx      int        `json:"nx"`      // X 方向的格点数
    Ny      int        `json:"ny"`      // Y 方向的格点数
}
```

#### DataStats
输入数据的统计信息。

```go
type DataStats struct {
    OriginalCount int     // 原始数据点数
    CleanedCount  int     // 清洗后有效数据点数
    NanCount      int     // NaN 值数量
    InfCount      int     // Inf 值数量
    DupCount      int     // 重复点数量
    MinX, MaxX    float64 // X 坐标范围
    MinY, MaxY    float64 // Y 坐标范围
    MinZ, MaxZ    float64 // Z 值范围
    MeanZ         float64 // Z 值均值
    StdZ          float64 // Z 值标准差
    IsCollinear   bool    // 数据是否近似共线
}
```

### 主要函数

#### GenerateLines — 等值线生成主入口
从散点数据生成等值线。完整流程：数据统计 → 参数校验 → 克里金插值 → 等值线提取。

```go
func GenerateLines(points []*Point, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature
```

#### GenerateLinesByGrid — 基于已有网格生成等值线
适用于需要自定义插值参数或复用已有网格的场景。

```go
func GenerateLinesByGrid(grid *Grid, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature
```

#### GenerateFaces — 等值面生成主入口
从等值线 features 生成等值面（填充多边形）。流程：解析等值线 → 递归切割多边形 → 去重 → 并行计算 level。

```go
func GenerateFaces(features []*geojson.Feature, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature
```

#### GenerateGrid — 克里金插值生成网格
执行 Ordinary Kriging 插值，失败时自动回退到 IDW。

```go
func GenerateGrid(points []*Point, boundary orb.Polygon, resolution int, levels []float64) (grid *Grid, err error)
```

#### ExtractPeaks — 提取山顶点
从插值网格中提取所有局部极大值点，按高度升序排列。

```go
func ExtractPeaks(grid *Grid, boundary orb.Polygon) []*Point
```

#### CombineGrid — 合并多个网格
通过自定义转换函数和聚合函数合并多个同尺寸网格，支持并行处理。

```go
func CombineGrid(gridList []*Grid, convFn func(gridIdx int, gridVal float64) float64, combineFn func(vals []float64) float64) (*Grid, error)
```

#### IdwFallBack — IDW 回退插值
当克里金插值失败时使用反距离加权法生成网格，结果归一化到原始 Z 值范围。

```go
func IdwFallBack(x, y, z, gridX, gridY []float64, levels []float64) (grid *Grid, err error)
```

### 校验函数

#### ValidatePoints — 数据点清洗和校验
过滤 nil、NaN、Inf、重复点，检测共线性。

```go
func ValidatePoints(points []*Point, minPoints int) (cleaned []*Point, err error)
```

#### ValidateBoundary — 边界多边形校验
检查多边形是否为空、顶点数、自交、面积。

```go
func ValidateBoundary(boundary orb.Polygon) error
```

#### ValidateContourOption — 配置参数校验
校验 Resolution 范围和 LevelList/ContourInterval 配置。

```go
func ValidateContourOption(opt *ContourOption) error
```

#### CheckParams — 完整参数校验
串联调用 ValidatePoints、ValidateContourOption、ValidateBoundary，并补全 LevelList。

```go
func CheckParams(points []*Point, boundary orb.Polygon, opt *ContourOption) (cleanPoints []*Point, err error)
```

#### ComputeDataStats — 计算统计信息
不修改原始数据，返回数据的统计概况。

```go
func ComputeDataStats(points []*Point) DataStats
```

### Grid 方法

| 方法 | 说明 |
|------|------|
| `Dims() (c, r int)` | 返回网格的列数和行数 |
| `X(c int) float64` | 返回第 c 列的 X 坐标 |
| `Y(r int) float64` | 返回第 r 行的 Y 坐标 |
| `Z(c, r int) float64` | 返回 (c, r) 位置的插值 |
| `Sample(lon, lat float64) float64` | 双线性插值采样 |
| `Normalize(...)` | 将网格值线性映射到目标范围 |
| `MarshalJSON() / UnmarshalJSON()` | JSON 序列化/反序列化 |

### 内部模块

#### internal/kriging — 克里金插值核心

| 函数/类型 | 说明 |
|-----------|------|
| `OrdinaryKriging` | 普通克里金插值器，保存数据、变差函数、模型统计 |
| `OKConfig` | 克里金配置（模型、参数、坐标类型、各向异性等） |
| `NewOrdinaryKriging(x, y, z, cfg)` | 创建并初始化克里金插值器 |
| `Execute(style, xpts, ypts, mask)` | 执行插值（grid/points/masked 模式） |
| `ExecuteGrid(xpts, ypts)` | 网格插值便捷方法 |
| `ComputeExperimentalVariogram(...)` | 计算实验变差函数（分箱统计） |
| `FitVariogramModel(...)` | 拟合变差函数参数（黄金分割 + 坐标下降） |
| `Krige(...)` | 单点克里金插值 |
| `FindStatistics(...)` | 留一交叉验证计算模型质量统计 |

#### internal/kriging — 变差函数模型

| 模型 | 公式 | 参数 |
|------|------|------|
| `linear` | γ(d) = slope·d + nugget | slope, nugget |
| `power` | γ(d) = scale·d^exp + nugget | scale, exponent, nugget |
| `gaussian` | γ(d) = psill·(1-exp(-d²/effR²)) + nugget | psill, range, nugget |
| `exponential` | γ(d) = psill·(1-exp(-d/effR)) + nugget | psill, range, nugget |
| `spherical` | γ(d) = psill·(1.5r-0.5r³) + nugget (r=d/range) | psill, range, nugget |
| `hole-effect` | γ(d) = psill·(1-(1-d/effR)·exp(-d/effR)) + nugget | psill, range, nugget |

#### internal/utils — 通用工具

| 函数 | 说明 |
|------|------|
| `CloseRing(ring)` | 闭合多边形环（首尾点不同时追加首点） |
| `EqPoint(a, b, tol)` | 判断两点在容差范围内是否相等 |
| `MakeLinespace(start, end, n)` | 生成 n 个均匀分布的浮点数 |
| `MinMax(items...)` | 返回切片的最小值和最大值（泛型） |
| `Min(items...)` | 返回切片的最小值（泛型） |
| `MatToFlat(a, rows, cols)` | 二维切片展平为一维（行优先） |

## 算法原理

### 克里金插值 (Ordinary Kriging)
1. **计算实验变差函数**：对数据点对的距离和半方差进行分箱统计
2. **拟合理论变差函数**：使用黄金分割搜索 + 坐标下降法 + SoftL1 鲁棒损失拟合模型参数
3. **构建克里金方程组**：(n+1)×(n+1) 矩阵，包含拉格朗日乘子约束
4. **求解线性方程组**：支持高斯消元和伪逆两种策略
5. **生成规则网格**：对每个网格点执行克里金插值

### IDW 回退
当克里金插值失败或结果包含 NaN/Inf 时自动启用：
- 使用 power=2 的反距离加权
- 坐标归一化平衡 X、Y 方向影响
- 结果归一化到原始 Z 值范围

### 等值线提取 (Marching Squares)
1. 使用 gonum/plot 的 `NewContour` 在自定义 canvas 上绘制等值线
2. 从 `vg.Path` 提取坐标并转换回地理空间
3. **线段合并**：空间哈希 + 并查集合并端点相接的线段
4. **边界裁剪**：Sutherland-Hodgman 风格算法裁剪到边界多边形内
5. **方向统一**：封闭线确保 CCW，开放线通过构造辅助多边形确定方向
6. **高值侧判断**：多点左法向采样 + 多数投票确定 high_value_side

### 等值面生成 (Polygon Splitting)
1. 从等值线 features 中分离边界和等值线
2. 解析 contourEntry（闭环按面积降序，开线在后）
3. **递归切割**：闭环形成外环+孔洞结构，开线将多边形一分为二
4. **去重**：通过规范化环的字符串 key 去重
5. **Level 计算**：并行遍历 contourEntry，通过共享边判断多边形所属 level
6. 转换为 GeoJSON Feature（规范化环方向、四舍五入坐标）

## 配置参数

### Resolution（分辨率）
- **范围建议**: 50-500
- **说明**: 指定较长边的格点数，较短边按比例自动计算
- **影响**: 分辨率越高，结果越精细，但计算时间越长

### ContourInterval（等值线间隔）
- **要求**: 必须为正数
- **说明**: 等值线之间的数值间隔
- **自动计算**: 当 ContourStart/End 为 0 时，自动从数据 Z 值范围推算

### LevelList（级别列表）
- **优先级**: 高于 ContourStart/ContourEnd/ContourInterval
- **格式**: 浮点数切片，如 `[100.0, 110.0, 120.0]`
- **用途**: 显式指定要生成的等值线级别

### 变差函数模型选择
- **spherical**（推荐）：最常用的地统计模型，适合大多数场景
- **exponential**：适合空间相关性随距离指数衰减的数据
- **gaussian**：适合空间相关性平滑衰减的数据
- **linear**：简单线性模型，参数最少
- **power**：幂函数模型，指数范围 0<exp<2

## 错误处理

### 数据校验
- 数据点数量不足（默认最少 3 个）
- NaN/Inf/重复点自动过滤并输出警告
- 共线数据检测（PCA 协方差矩阵行列式）

### 边界校验
- 多边形为空或顶点数不足
- 外环自交检测（非相邻线段对严格相交判断）
- 面积过小检测

### 插值容错
- 克里金插值失败 → 自动回退到 IDW
- 矩阵奇异 → 高斯消元失败 → 回退到伪逆
- 网格值异常（NaN/Inf）→ 回退到 IDW

## 示例输出

### 等值线 GeoJSON

```json
{
  "type": "Feature",
  "geometry": {
    "type": "LineString",
    "coordinates": [[120.0, 30.0], [120.1, 30.1], ...]
  },
  "properties": {
    "type": "temperature",
    "level": 100.0,
    "is_closed": false,
    "high_value_side": "left",
    "source": "weather_station"
  }
}
```

### 等值面 GeoJSON

```json
{
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates": [[[120.0, 30.0], [120.1, 30.0], [120.1, 30.1], [120.0, 30.0]]]
  },
  "properties": {
    "type": "temperature",
    "level": 100.0,
    "source": "weather_station"
  }
}
```

### 图像输出
支持生成 1600×1200 像素的 PNG 图片：
- 热力图底色（按 level 分带着色，蓝→黄→橙→红渐变）
- 等值线叠加
- 坐标系为 (lat, lon)

## 依赖项

| 依赖 | 用途 |
|------|------|
| `github.com/panjf2000/ants/v2` | goroutine 池管理，并行处理 |
| `github.com/paulmach/orb` | 几何数据处理（Point、Ring、Polygon、GeoJSON） |
| `gonum.org/v1/gonum` | 矩阵运算（mat.Dense）、线性方程组求解 |
| `gonum.org/v1/plot` | 等值线绘制（Marching Squares）、图像渲染 |
| `github.com/duke-git/lancet/v2` | 工具函数（maputil.Merge） |

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request。

## 作者

[inkbamboo](https://github.com/inkbamboo)

## 版本历史

- v0.0.1: 初始版本，支持等值线生成
- v0.0.2: 添加等值面生成功能
- v0.0.3: 添加地形特征提取和图像渲染
