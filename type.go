package kriging_contour

// ContourOption 定义等值线/等值面的生成参数配置。
// 可通过 JSON 反序列化或直接在代码中构造。
type ContourOption struct {
	// ContourType 等值线类型名称，会作为属性写入输出的 GeoJSON Feature 中（如 "comprehensive"）
	ContourType string `json:"contour_type"`

	// ContourStart 等值线的起始数值。当 LevelList 为空时，与 ContourEnd、ContourInterval 配合生成 LevelList。
	ContourStart float64 `json:"contour_start"`

	// ContourEnd 等值线的结束数值。与 ContourStart、ContourInterval 配合生成 LevelList。
	ContourEnd float64 `json:"contour_end,omitempty"`

	// ContourInterval 等值线之间的数值间隔。用于从 Start/End 自动生成 LevelList。
	ContourInterval float64 `json:"contour_interval"`

	// LevelList 显式指定的等值线级别列表。若不为空，则优先使用，忽略 Start/End/Interval。
	LevelList []float64 `json:"level_list,omitempty"`

	// Resolution 克里金插值网格的分辨率（长边方向的格点数）。范围建议 50-500。
	Resolution int `json:"resolution"`

	// ImagePath 可选输出等值线渲染图片的路径。为空则不输出图片。
	ImagePath string `json:"image_path,omitempty"`

	// Extra 附加属性，将合并到每个输出的 GeoJSON Feature 的 properties 中。
	Extra map[string]interface{} `json:"extra"`
}

// Point 表示一个空间采样数据点。
// X 为经度或 X 坐标，Y 为纬度或 Y 坐标，Z 为该点的观测值。
type Point struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
	Z float64 `json:"z"`
}
