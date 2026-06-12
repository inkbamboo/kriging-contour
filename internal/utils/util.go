// Package utils 提供 GeoJSON 读写、几何计算、数学工具等通用辅助函数。
package utils

import (
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/bytedance/sonic"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	"github.com/paulmach/orb/planar"
	"golang.org/x/exp/constraints"
)

// LoadGeoJSON 从指定文件路径加载GeoJSON数据并解析为要素集合
// 参数:
//   - geoFile: GeoJSON文件路径
//
// 返回值:
//   - *geojson.FeatureCollection: 解析后的要素集合
//   - error: 读取或解析失败时返回错误
func LoadGeoJSON(geoFile string) (geoData *geojson.FeatureCollection, err error) {
	var data []byte
	if data, err = os.ReadFile(geoFile); err != nil {
		err = fmt.Errorf("读取Geo文件失败: %v\n", err)
		return
	}
	geoData, err = geojson.UnmarshalFeatureCollection(data)
	if err != nil {
		return nil, fmt.Errorf("解析边界GeoJSON失败: %w", err)
	}
	return geoData, nil
}

// WriteGeoJson 将GeoJSON要素集合序列化后写入指定文件，自动创建父目录
// 参数:
//   - geoData: 待写入的要素集合
//   - filePath: 输出文件路径
//
// 返回值:
//   - error: 序列化或写入失败时返回错误
func WriteGeoJson(geoData *geojson.FeatureCollection, filePath string) (err error) {
	fileDir := filepath.Dir(filePath)
	if fileDir != "." {
		if err = os.MkdirAll(fileDir, 0755); err != nil {
			return
		}
	}
	var data []byte
	if data, err = sonic.Marshal(geoData); err != nil {
		return err
	}
	return os.WriteFile(filePath, data, 0644)
}

// Max 使用泛型计算输入参数中的最大值
// 参数:
//   - items: 可比较类型的参数列表
//
// 返回值:
//   - T: 最大值
func Max[T constraints.Ordered](items ...T) T {
	if len(items) == 0 {
		return *new(T)
	}
	maxVal := items[0]
	for _, v := range items {
		if maxVal < v {
			maxVal = v
		}
	}
	return maxVal
}

// Distance 计算两个坐标点之间的欧几里得距离
// 参数:
//   - a, b: 两个坐标点
//
// 返回值:
//   - float64: 欧几里得距离
func Distance(a, b orb.Point) float64 {
	dx := a[0] - b[0]
	dy := a[1] - b[1]
	return math.Sqrt(dx*dx + dy*dy)
}

// ReverseCoords 反转坐标点列表的顺序（不修改原切片）
// 参数:
//   - coords: 原始坐标点列表
//
// 返回值:
//   - []orb.Point: 反转后的新坐标点列表
func ReverseCoords(coords []orb.Point) []orb.Point {
	n := len(coords)
	result := make([]orb.Point, n)
	for i := 0; i < n; i++ {
		result[i] = coords[n-1-i]
	}
	return result
}

// PolyEqual 判断两个多边形是否近似相等（通过环长度和面积差比较）
// 参数:
//   - a, b: 待比较的两个多边形
//
// 返回值:
//   - bool: 面积差小于1.0时视为相等
func PolyEqual(a, b orb.Polygon) bool {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		return false
	}
	areaDiff := math.Abs(planar.Area(a) - planar.Area(b))
	return areaDiff < 1.0
}

// PolyArea 计算多边形的面积（便捷封装 planar.Area）
// 参数:
//   - poly: 多边形
//
// 返回值:
//   - float64: 多边形面积
func PolyArea(poly orb.Polygon) float64 {
	return planar.Area(poly)
}

// IsClosedLine 判断坐标点列表是否构成闭合线（首尾距离小于容差）
// 参数:
//   - coords: 坐标点列表
//   - tol: 闭合判定容差
//
// 返回值:
//   - bool: 首尾距离小于tol则为闭合线
func IsClosedLine(coords []orb.Point, tol float64) bool {
	if len(coords) < 3 {
		return false
	}
	return Distance(coords[0], coords[len(coords)-1]) < tol
}

// CloseRing 确保多边形环的首尾点相同（闭合），如果不同则追加首点到末尾
// 参数:
//   - ring: 多边形环的坐标点列表
//
// 返回值:
//   - orb.Ring: 首尾闭合的多边形环
func CloseRing(ring orb.Ring) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	if ring[0] == ring[len(ring)-1] {
		return ring
	}
	closed := make(orb.Ring, len(ring)+1)
	copy(closed, ring)
	closed[len(ring)] = ring[0]
	return closed
}

// Round6 将浮点数四舍五入保留6位小数
// 参数:
//   - v: 原始浮点数
//
// 返回值:
//   - float64: 保留6位小数的结果
func Round6(v float64) float64 {
	return math.Round(v*1e6) / 1e6
}

// ExtractArc 从多边形环上按起止距离提取弧段
// 支持正常弧段（startD ≤ endD）和跨越末端的弧段（startD > endD）
// 参数:
//   - ring: 多边形环的坐标点列表
//   - startD: 弧段起点在环上的累计距离
//   - endD: 弧段终点在环上的累计距离
//
// 返回值:
//   - orb.LineString: 提取的弧段坐标点列表
func ExtractArc(ring orb.Ring, startD, endD float64) orb.LineString {
	var points []orb.Point
	cumDist := 0.0

	collecting := false
	if startD <= endD {
		// 正常弧段
		for i := 0; i < len(ring); i++ {
			if i == 0 {
				if cumDist <= startD && startD < cumDist+Distance(ring[0], ring[1]) {
					// 插入start点
					t := (startD - cumDist) / Distance(ring[0], ring[1])
					p := orb.Point{
						ring[0][0] + t*(ring[1][0]-ring[0][0]),
						ring[0][1] + t*(ring[1][1]-ring[0][1]),
					}
					points = append(points, p)
					collecting = true
				}
			}
			if collecting {
				points = append(points, ring[i])
			}
			if i < len(ring)-1 {
				segLen := Distance(ring[i], ring[i+1])
				if collecting && cumDist+segLen >= endD {
					t := (endD - cumDist) / segLen
					p := orb.Point{
						ring[i][0] + t*(ring[i+1][0]-ring[i][0]),
						ring[i][1] + t*(ring[i+1][1]-ring[i][1]),
					}
					points = append(points, p)
					break
				}
				cumDist += segLen
			}
		}
	} else {
		// 跨越末端: startD -> end, 0 -> endD
		// 先 startD -> blen
		collecting = false
		for i := 0; i < len(ring); i++ {
			if i == 0 && cumDist <= startD && startD < cumDist+Distance(ring[0], ring[1]) {
				t := (startD - cumDist) / Distance(ring[0], ring[1])
				p := orb.Point{
					ring[0][0] + t*(ring[1][0]-ring[0][0]),
					ring[0][1] + t*(ring[1][1]-ring[0][1]),
				}
				points = append(points, p)
				collecting = true
			}
			if collecting {
				points = append(points, ring[i])
			}
			if i < len(ring)-1 {
				segLen := Distance(ring[i], ring[i+1])
				cumDist += segLen
			}
		}
		// 然后 0 -> endD
		cumDist = 0
		for i := 0; i < len(ring); i++ {
			if cumDist+Distance(ring[i], ring[i+1]) >= endD {
				t := (endD - cumDist) / Distance(ring[i], ring[i+1])
				p := orb.Point{
					ring[i][0] + t*(ring[i+1][0]-ring[i][0]),
					ring[i][1] + t*(ring[i+1][1]-ring[i][1]),
				}
				points = append(points, p)
				break
			}
			points = append(points, ring[i])
			cumDist += Distance(ring[i], ring[i+1])
		}
	}

	if len(points) < 2 {
		return orb.LineString{ring[0], ring[0]}
	}
	return points
}

// ProjectToRing 计算点在多边形环上的最近投影距离（沿环的累计距离）
// 参数:
//   - pt: 待投影的坐标点
//   - ring: 多边形环的坐标点列表
//
// 返回值:
//   - float64: 投影点在环上的累计距离
func ProjectToRing(pt orb.Point, ring orb.Ring) float64 {
	minDist := math.Inf(1)
	bestD := 0.0
	cumDist := 0.0

	for i := 0; i < len(ring)-1; i++ {
		a, b := ring[i], ring[i+1]
		// 计算点到线段投影
		dx := b[0] - a[0]
		dy := b[1] - a[1]
		lenSq := dx*dx + dy*dy

		var proj orb.Point
		if lenSq < 1e-12 {
			proj = a
		} else {
			t := ((pt[0]-a[0])*dx + (pt[1]-a[1])*dy) / lenSq
			if t < 0 {
				proj = a
			} else if t > 1 {
				proj = b
			} else {
				proj = orb.Point{a[0] + t*dx, a[1] + t*dy}
			}
		}

		d := Distance(pt, proj)
		if d < minDist {
			minDist = d
			segDist := math.Sqrt(dx*dx + dy*dy)
			if lenSq > 1e-12 {
				t := ((pt[0]-a[0])*dx + (pt[1]-a[1])*dy) / lenSq
				if t < 0 {
					t = 0
				} else if t > 1 {
					t = 1
				}
				bestD = cumDist + t*segDist
			} else {
				bestD = cumDist
			}
		}
		cumDist += math.Sqrt(dx*dx + dy*dy)
	}
	return bestD
}

// MakeLinspace 生成从 start 到 end 的 n 个等间距点。
// 若 n <= 1，返回 [start]。
func MakeLinspace(start, end float64, n int) []float64 {
	if n <= 1 {
		return []float64{start}
	}
	result := make([]float64, n)
	step := (end - start) / float64(n-1)
	for i := 0; i < n; i++ {
		result[i] = start + float64(i)*step
	}
	return result
}

// MinMax 使用泛型计算输入参数的最小值和最大值
// 参数:
//   - items: 可比较类型的参数列表
//
// 返回值:
//   - T: 最小值
//   - T: 最大值
func MinMax[T constraints.Ordered](items ...T) (T, T) {
	if len(items) == 0 {
		return *new(T), *new(T)
	}
	minVal, maxVal := items[0], items[0]
	for _, v := range items[1:] {
		if v < minVal {
			minVal = v
		}
		if v > maxVal {
			maxVal = v
		}
	}
	return minVal, maxVal
}
