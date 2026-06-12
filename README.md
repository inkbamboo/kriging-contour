# kriging-contour

基于 [克里金插值（Ordinary Kriging）](https://en.wikipedia.org/wiki/Kriging) 的等值线/等值面生成库，Go 语言实现。

## 概述

从离散空间数据点出发，通过普通克里金插值生成规则网格，再使用 Marching Squares 算法提取等值线（LineString），最后通过遍历分割生成等值面（Polygon），所有结果输出为标准 GeoJSON 格式。

```
离散数据点 → 克里金插值 → 规则网格 → Marching Squares → 等值线 → 遍历分割 → 等值面
```

## 安装

```bash
go get github.com/inkbamboo/kriging-contour
```

## 快速开始

```go
package main

import (
    kc "github.com/inkbamboo/kriging-contour"
    "github.com/inkbamboo/kriging-contour/internal/utils"
    "github.com/paulmach/orb"
    "github.com/paulmach/orb/geojson"
)

func main() {
    // 1. 加载边界多边形
    boundaryFC, _ := utils.LoadGeoJSON("boundary.json")
    boundary := boundaryFC.Features[0].Geometry.(orb.Polygon)

    // 2. 准备数据点 (X=lon, Y=lat, Z=测量值)
    points := []*kc.Point{
        {X: 116.0, Y: 39.0, Z: 25.5},
        {X: 116.5, Y: 39.0, Z: 30.2},
    }

    // 3. 配置参数
    opt := kc.ContourOption{
        ContourType:     "temperature",
        ContourStart:    20,
        ContourEnd:      40,
        ContourInterval: 5,
        Resolution:      200,
        ImagePath:       "contour.png",
        Extra: map[string]interface{}{"unit": "℃"},
    }

    // 4. 生成等值线
    lineFeatures := kc.GenerateLines(points, boundary, opt)

    // 5. 生成等值面
    faceFeatures := kc.GenerateFaces(lineFeatures, boundary, opt)

    // 6. 导出 GeoJSON
    lineFC := geojson.NewFeatureCollection()
    lineFC.Features = lineFeatures
    utils.WriteGeoJson(lineFC, "output_lines.json")

    faceFC := geojson.NewFeatureCollection()
    faceFC.Features = faceFeatures
    utils.WriteGeoJson(faceFC, "output_faces.json")
}
```

## API 参考

### 数据类型

| 类型 | 说明 |
|------|------|
| `Point` | 三维数据点，X=经度, Y=纬度, Z=测量值 |
| `GridData` | 克里金插值网格，Z 为 ny×nx 矩阵 |
| `ContourOption` | 等值线/面生成参数配置 |

### ContourOption 参数

| 字段 | 类型 | 必填 | 说明 |
|------|------|------|------|
| `ContourType` | string | 否 | 类型标识，透传到 GeoJSON properties["type"] |
| `ContourStart` | float64 | 是 | 等值线起始值 |
| `ContourEnd` | float64 | 否 | 等值线结束值，0 表示自动推断 |
| `ContourInterval` | float64 | 是 | 等值线级别间隔 |
| `Resolution` | int | 是 | 克里金网格分辨率（长边方向网格数） |
| `ImagePath` | string | 否 | 等值线热力图 PNG 输出路径 |
| `Extra` | map | 否 | 透传到 GeoJSON Properties 的额外字段 |

### 核心函数

| 函数 | 说明 |
|------|------|
| `GenerateLines(points, boundary, opt) []*geojson.Feature` | 生成等值线（LineString），含裁剪、合并、方向统一、短线过滤 |
| `GenerateFaces(features, boundary, opt) []*geojson.Feature` | 基于等值线生成等值面（Polygon） |

### 等值线 Feature Properties

| 字段 | 类型 | 说明 |
|------|------|------|
| `level` | float64 | 等值线级别值 |
| `is_closed` | bool | 是否闭合 |
| `high_value_side` | string | 高值侧：`left`/`right`（开线）, `inside`/`outside`（闭合环） |
| `type` | string | `"boundary"` 表示边界，否则为透传的 ContourType |

### 等值面 Feature Properties

| 字段 | 类型 | 说明 |
|------|------|------|
| `level` | float64 | 等值面级别值 |
| `type` | string | 透传的 ContourType |

## 内部实现

### GenerateLines — 等值线生成

1. **克里金插值** (`generateGrid`)：用球状变异函数模型的普通克里金生成规则网格并归一化缩放
2. **级别计算** (`getContourLevels`)：按 ContourStart/End/Interval 生成 level 列表
3. **并发提取**：用 goroutine 池为每个 level 并行执行 marching squares，通过记录型 Canvas 截获路径
4. **坐标映射** (`extractCoords`)：画布坐标逆映射到数据坐标系
5. **边界裁剪** (`clipLineToPolygon`)：保留边界内部部分，处理四种出入情况
6. **线段合并** (`mergeConnectedLines`)：四种连接方式迭代合并直到收敛
7. **短线过滤**：丢弃长度 < 边界周长 × 0.01 的线
8. **方向统一** (`resolveLineOrientationCG`)：右手定则 + 坐标转换后 CCW/CW 翻转

### GenerateFaces — 等值面生成

1. **解析** (`parseContourEntries`)：提取闭合环和开线，闭环按面积从大到小排序
2. **遍历分割** (`splitPolys`)：对每条等值线依次切分所有已产生的子多边形
3. **去重** (`deduplicatePolygons`)：通过归一化环的序列化键去重
4. **并发计算 Level** (`computePolygonLevel`)：goroutine 池并发为每个面推断 level 值

### 关键几何算法

| 函数 | 说明 |
|------|------|
| `segmentIntersection` | 参数方程法求两线段交点 |
| `clipLineToPolygon` | Sutherland-Hodgman 风格线裁剪 |
| `splitRingAtPoints` | 环在两点处拆分为两段弧 |
| `resolveLineOrientationCG` | 右手定则方向统一（双线性采样判断高值侧） |
| `signedRingArea` | 鞋带公式计算有符号面积 |

## 变异函数模型

克里金插值支持以下变异函数模型，参数可自动拟合：

| 模型 | 参数 |
|------|------|
| `linear` | [slope, nugget] |
| `power` | [scale, exponent, nugget] |
| `gaussian` | [psill, range, nugget] |
| `spherical` | [psill, range, nugget]（默认） |
| `exponential` | [psill, range, nugget] |
| `hole-effect` | [psill, range, nugget] |

## 坐标系统与方向约定

- 输入 Point：`[lon, lat]`
- 边界多边形：`[lat, lon]`
- 内部计算：`[lon, lat]`
- GeoJSON 输出：`[lat, lon]`

等值线遵循**右手定则**：沿线条走向，高值区在左侧。闭合环按逆时针（CCW）。

## 项目结构

```
kriging-contour/
├── contour.go              # 网格生成、等值线提取、裁剪、合并、热力图
├── faces.go                # 等值面生成（解析→遍历分割→去重→Level计算）
├── type.go                 # ContourOption / Point / GridData
├── internal/
│   ├── contour/            # Grid 适配器 + Canvas 记录器
│   ├── kriging/            # 克里金插值（PyKrige 移植）
│   │   ├── ok.go           # OrdinaryKriging 类型
│   │   ├── core.go         # 距离计算、变异函数拟合、求解
│   │   ├── grid.go         # 网格坐标生成
│   │   └── variogram.go    # 6 种变异函数模型
│   └── utils/              # GeoJSON 读写、几何计算
└── files/                  # 测试数据
```

## 依赖

| 依赖 | 用途 |
|------|------|
| [gonum](https://github.com/gonum/gonum) | 矩阵运算 |
| [gonum/plot](https://github.com/gonum/plot) | Marching Squares + 热力图 |
| [paulmach/orb](https://github.com/paulmach/orb) | GIS 几何类型 + GeoJSON |
| [panjf2000/ants](https://github.com/panjf2000/ants) | Goroutine 池 |
| [lancet](https://github.com/duke-git/lancet) | 工具函数 |
| [sonic](https://github.com/bytedance/sonic) | 高性能 JSON |

## 参考

- [PyKrige](https://github.com/GeoStat-Framework/PyKrige) — Python 克里金插值库
- P.K. Kitanidis, *Introduction to Geostatistics*, Cambridge, 1997

## License

MIT
