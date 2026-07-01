# Kriging Contour - 克里金等值线/等值面生成库

一个基于克里金插值（Kriging）的等值线（Contour）和等值面（Contour Faces）生成库，支持从离散采样点数据生成高质量的 GeoJSON 格式等值线和等值面。

## 功能特性

### 核心功能
- **等值线生成**：从离散点数据生成等值线（Contour Lines）
- **等值面生成**：从等值线生成填充多边形等值面（Contour Faces）
- **克里金插值**：使用 Ordinary Kriging 进行空间插值
- **IDW 回退**：当克里金插值失败时自动回退到反距离加权插值
- **地形特征提取**：自动检测山顶/洼地等特征点

### 数据处理
- **数据清洗**：自动过滤 NaN、Inf 和重复数据点
- **参数验证**：完整的输入参数和边界校验
- **统计信息**：计算数据范围、均值、标准差等统计量

### 输出格式
- **GeoJSON 输出**：标准 GeoJSON FeatureCollection 格式
- **图像渲染**：可选的等值线热力图 PNG 输出
- **网格序列化**：支持网格数据的 JSON 序列化/反序列化

### 性能优化
- **并行处理**：使用 goroutine 池（ants）并行处理
- **内存优化**：预分配切片，避免重复内存分配
- **缓存优化**：行优先遍历，提高缓存命中率

## 项目结构

```
kriging-contour/
├── type.go              # 数据结构和类型定义
├── contour.go           # 等值线生成主逻辑
├── faces.go             # 等值面生成算法
├── grid.go              # 网格数据结构和克里金插值
├── validate.go          # 数据校验和参数验证
├── peaks.go             # 地形特征提取
├── image.go             # 等值线图像渲染
├── canvas.go            # 自定义 Canvas 实现
├── go.mod              # Go 模块依赖
└── internal/           # 内部工具库
    ├── kriging/        # 克里金插值实现
    └── utils/          # 通用工具函数
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
    ContourType      string                 `json:"contour_type"`      // 等值线类型名称
    ContourStart     float64                `json:"contour_start"`     // 起始数值
    ContourEnd       float64                `json:"contour_end,omitempty"` // 结束数值
    ContourInterval  float64                `json:"contour_interval"`  // 数值间隔
    LevelList        []float64              `json:"level_list,omitempty"` // 显式级别列表
    Resolution       int                    `json:"resolution"`        // 网格分辨率
    ImagePath        string                 `json:"image_path,omitempty"` // 图片输出路径
    Extra            map[string]interface{} `json:"extra"`             // 附加属性
}
```

#### Grid
克里金插值后的规则网格数据结构。

```go
type Grid struct {
    Data    *mat.Dense // 插值后的网格数据矩阵
    XCoords []float64  `json:"xCoords"` // X 坐标轴（列坐标）
    YCoords []float64  `json:"yCoords"` // Y 坐标轴（行坐标）
    Nx      int        `json:"nx"`      // X 方向的格点数
    Ny      int        `json:"ny"`      // Y 方向的格点数
}
```

### 主要函数

#### GenerateLines
从散点数据生成等值线的主入口函数。

```go
func GenerateLines(points []*Point, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature
```

**参数**:
- `points`: 离散采样点数据
- `boundary`: 边界多边形，用于裁剪等值线和限制插值范围
- `opt`: 等值线生成配置参数

**返回**: GeoJSON Feature 切片，每个 Feature 代表一条等值线

#### GenerateFaces
从等值线 features 生成等值面（填充多边形）。

```go
func GenerateFaces(features []*geojson.Feature, boundary orb.Polygon, opt *ContourOption) []*geojson.Feature
```

**参数**:
- `features`: `GenerateLines` 返回的等值线 GeoJSON Feature 列表（含边界 Feature）
- `boundary`: 原始边界多边形
- `opt`: 等值线配置参数

**返回**: 等值面 GeoJSON Feature 切片

#### GenerateGrid
执行克里金插值生成规则网格。

```go
func GenerateGrid(points []*Point, boundary orb.Polygon, resolution int, levels []float64) (grid *Grid, err error)
```

**参数**:
- `points`: 离散采样点数据
- `boundary`: 边界多边形
- `resolution`: 网格分辨率（长边方向的格点数）
- `levels`: 等值线级别列表

**返回**: 插值后的网格数据和可能的错误

#### ExtractPeaks
从插值网格中提取所有山顶的位置及高度。

```go
func ExtractPeaks(grid *Grid, boundary orb.Polygon) []*Point
```

**参数**:
- `grid`: 插值网格
- `boundary`: 边界多边形

**返回**: 山顶点切片，按高度升序排列

#### CombineGrid
将多个网格数据合并为单个网格。

```go
func CombineGrid(gridList []*Grid, convFn func(gridIdx int, gridVal float64) float64, combineFn func(vals []float64) float64) (*Grid, error)
```

**参数**:
- `gridList`: 待合并的网格切片
- `convFn`: 转换函数，接收网格索引和该网格在当前位置的值，返回转换后的值
- `combineFn`: 聚合函数，接收每个网格转换后的值切片，返回合并结果

**返回**: 合并后的网格和可能的错误

## 算法原理

### 克里金插值
使用 Ordinary Kriging 进行空间插值：
1. 计算实验变异函数
2. 拟合理论变异函数模型
3. 求解克里金方程组
4. 生成规则网格插值结果

### IDW 回退
当克里金插值失败或结果异常时，自动回退到反距离加权插值：
- 使用 power=2 的权重函数
- 坐标归一化以平衡 X、Y 方向的影响
- 结果归一化到原始数据的 Z 值范围

### 等值线提取
基于 Marching Squares 算法：
1. 使用 gonum/plot 的 NewContour 在自定义 canvas 上绘制
2. 从 vg.Path 提取坐标并转换回地理空间
3. 合并端点相接的线段
4. 用边界多边形裁剪
5. 统一线段方向（resolveLineOrientationCG）

### 等值面生成
基于多边形切割算法：
1. 从等值线中分离边界特征和等值线特征
2. 递归用每条等值线切割边界多边形
3. 去重（deduplicatePolygons）
4. 并行计算每个面的 level
5. 转换为 GeoJSON Feature

## 性能优化

### 并行处理
- 使用 `github.com/panjf2000/ants/v2` 进行 goroutine 池管理
- 并行处理每个等值线级别的提取
- 并行计算每个等值面的 level

### 内存优化
- 预分配切片，避免重复内存分配
- 使用引用传递，避免不必要的数据拷贝
- 行优先遍历，提高缓存命中率

### 算法优化
- 空间哈希加速邻近搜索
- 并查集（Union-Find）合并连通分量
- 二分查找定位网格单元

## 配置参数

### Resolution（分辨率）
- **范围建议**: 50-500
- **说明**: 指定较长边的格点数
- **影响**: 分辨率越高，结果越精细，但计算时间越长

### ContourInterval（等值线间隔）
- **要求**: 必须为正数
- **说明**: 等值线之间的数值间隔
- **示例**: 10.0 表示每隔 10 个单位绘制一条等值线

### LevelList（级别列表）
- **优先级**: 高于 ContourStart/ContourEnd/ContourInterval
- **格式**: 浮点数切片，如 `[100.0, 110.0, 120.0]`
- **用途**: 显式指定要生成的等值线级别

## 错误处理

库提供了完整的错误处理机制：

### 数据校验错误
- 数据点数量不足
- 边界多边形无效（自交、面积过小）
- 参数配置错误

### 插值错误
- 克里金插值失败（自动回退到 IDW）
- 网格尺寸不一致
- 内存分配失败

### 几何处理错误
- 等值线提取失败
- 多边形切割失败
- 坐标转换错误

## 示例输出

### GeoJSON 格式

```json
{
  "type": "FeatureCollection",
  "features": [
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
    },
    // ... 更多等值线
  ]
}
```

### 图像输出
支持生成 1600×1200 像素的 PNG 图片，包含：
- 热力图底色（按 level 分带着色）
- 等值线叠加
- 从冷色（蓝）到暖色（红）的渐变色带

## 依赖项

- `github.com/panjf2000/ants/v2`: goroutine 池管理
- `github.com/paulmach/orb`: 几何数据处理
- `gonum.org/v1/gonum`: 矩阵运算
- `gonum.org/v1/plot`: 等值线绘制和图像渲染
- `github.com/duke-git/lancet/v2`: 工具函数

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