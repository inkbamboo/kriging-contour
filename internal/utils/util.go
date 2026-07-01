package utils

import (
	"math"

	"github.com/paulmach/orb"
	"golang.org/x/exp/constraints"
)

// CloseRing 闭合一个环（Ring）。如果首尾点不相等，则将首点追加到末尾。
// 用于确保 GeoJSON 多边形环符合规范要求（首尾闭合）。
func CloseRing(ring orb.Ring) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	if !ring[0].Equal(ring[len(ring)-1]) {
		ring = append(ring, ring[0])
	}
	return ring
}

// EqPoint 判断两点是否在给定的容差范围内相等。
// 使用欧几里得距离进行计算。
func EqPoint(a, b orb.Point, tolerance float64) bool {
	return math.Hypot(a[0]-b[0], a[1]-b[1]) <= tolerance
}

// MakeLinespace 在 [start, end] 区间内生成 n 个均匀分布的浮点数。
// 如果 n <= 1，则返回只包含 start 的切片。
func MakeLinespace(start, end float64, n int) []float64 {
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

// MinMax 返回切片元素中的最小值和最大值。
// 使用泛型支持任何 constraints.Ordered 类型（如 float64、int 等）。
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

func Min[T constraints.Ordered](items ...T) T {
	if len(items) == 0 {
		return *new(T)
	}
	minVal := items[0]
	for _, v := range items[1:] {
		if v < minVal {
			minVal = v
		}
	}
	return minVal
}

// MatToFlat 将二维浮点数切片展平为一维切片（行优先顺序）。
// 用于将 [][]float64 矩阵转换为 gonum/mat.Dense 所需的一维数据格式。
func MatToFlat(a [][]float64, rows, cols int) []float64 {
	flat := make([]float64, rows*cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			flat[i*cols+j] = a[i][j]
		}
	}
	return flat
}
