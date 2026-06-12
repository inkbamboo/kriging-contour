package kriging_contour

import "gonum.org/v1/gonum/mat"

// ContourOption 等值线/等值面生成参数配置
type ContourOption struct {
	// ContourType 等值线类型标识，透传到输出 GeoJSON properties["type"]
	ContourType string `json:"contour_type"`
	// ContourStart 等值线起始值（最小值）
	ContourStart float64 `json:"contour_start"`
	// ContourEnd 等值线结束值（最大值），0 表示根据数据自动推断
	ContourEnd float64 `json:"contour_end,omitempty"`
	// ContourInterval 等值线级别间隔
	ContourInterval float64 `json:"contour_interval"`
	// Resolution 克里金插值网格分辨率（长边方向的网格数）
	Resolution int `json:"resolution"`
	// ImagePath 可选：等值线热力图 PNG 输出路径
	ImagePath string `json:"image_path,omitempty"`
	// Extra 额外参数，会透传到输出的 GeoJSON Properties 中
	Extra map[string]interface{} `json:"extra"`
}

// Point 表示一个三维空间数据点，用于克里金插值输入。
// X/Y 为经纬度坐标（[lon, lat]），Z 为对应的测量值。
type Point struct {
	X float64 `json:"x"` // 经度（longitude）
	Y float64 `json:"y"` // 纬度（latitude）
	Z float64 `json:"z"` // 测量值
}

// GridData 表示克里金插值后的网格数据。
// Z 矩阵按行优先存储，行对应 Y（纬度），列对应 X（经度）。
type GridData struct {
	Z *mat.Dense // 插值结果矩阵，ny 行 × nx 列
	X []float64  // X 轴（经度）坐标序列，长度 = nx
	Y []float64  // Y 轴（纬度）坐标序列，长度 = ny
}
