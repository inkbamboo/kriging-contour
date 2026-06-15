// Package utils 日志工具，提供分步骤的结构化日志输出。
package utils

import (
	"fmt"
	"log"
	"os"
	"time"
)

// LogLevel 日志级别
type LogLevel int

const (
	LogLevelDebug LogLevel = iota
	LogLevelInfo
	LogLevelWarn
	LogLevelError
)

// StepLogger 步骤日志记录器，记录每个处理步骤的起止时间、状态和错误。
type StepLogger struct {
	logger *log.Logger
	level  LogLevel
	steps  []StepRecord
}

// StepRecord 单步执行记录
type StepRecord struct {
	Name     string
	Start    time.Time
	Duration time.Duration
	Status   string // "ok", "warn", "error", "skipped"
	Err      error
	Extra    map[string]interface{}
}

// NewStepLogger 创建一个新的步骤日志记录器。
// output: 日志输出目标，nil 表示 os.Stdout
func NewStepLogger(output *os.File) *StepLogger {
	if output == nil {
		output = os.Stdout
	}
	return &StepLogger{
		logger: log.New(output, "[kriging-contour] ", log.Ltime|log.Lmicroseconds),
		level:  LogLevelInfo,
	}
}

// SetLevel 设置日志级别
func (l *StepLogger) SetLevel(level LogLevel) {
	l.level = level
}

// Debug 输出 DEBUG 级别日志
func (l *StepLogger) Debug(format string, args ...interface{}) {
	if l.level <= LogLevelDebug {
		l.logger.Printf("[DEBUG] "+format, args...)
	}
}

// Info 输出 INFO 级别日志
func (l *StepLogger) Info(format string, args ...interface{}) {
	if l.level <= LogLevelInfo {
		l.logger.Printf("[INFO]  "+format, args...)
	}
}

// Warn 输出 WARN 级别日志
func (l *StepLogger) Warn(format string, args ...interface{}) {
	if l.level <= LogLevelWarn {
		l.logger.Printf("[WARN]  "+format, args...)
	}
}

// Error 输出 ERROR 级别日志
func (l *StepLogger) Error(format string, args ...interface{}) {
	if l.level <= LogLevelError {
		l.logger.Printf("[ERROR] "+format, args...)
	}
}

// BeginStep 开始一个新步骤并返回一个可 defer 的函数来自动结束步骤。
// 用法:
//
//	end := logger.BeginStep("数据校验")
//	defer end(nil, nil)
//	// ... 执行步骤 ...
//	end(nil, extra)  // 或 end(err, extra)
func (l *StepLogger) BeginStep(name string) func(err error, extra map[string]interface{}) StepRecord {
	rec := StepRecord{
		Name:  name,
		Start: time.Now(),
	}
	l.Info("步骤开始: %s", name)
	l.steps = append(l.steps, rec)

	return func(err error, extra map[string]interface{}) StepRecord {
		rec.Duration = time.Since(rec.Start)
		rec.Err = err
		rec.Extra = extra

		switch {
		case err != nil:
			rec.Status = "error"
			l.Error("步骤失败: %s (耗时: %v) - %v", name, rec.Duration.Round(time.Millisecond), err)
		case extra != nil && extra["warning"] != nil:
			rec.Status = "warn"
			l.Warn("步骤完成(有警告): %s (耗时: %v)", name, rec.Duration.Round(time.Millisecond))
		default:
			rec.Status = "ok"
			l.Info("步骤完成: %s (耗时: %v)", name, rec.Duration.Round(time.Millisecond))
		}

		// 更新记录
		for i := range l.steps {
			if l.steps[i].Name == name && l.steps[i].Start.Equal(rec.Start) {
				l.steps[i] = rec
				break
			}
		}
		return rec
	}
}

// Summary 输出所有步骤摘要
func (l *StepLogger) Summary() {
	fmt.Fprintln(l.logger.Writer(), "\n========== 处理摘要 ==========")
	errCount := 0
	warnCount := 0
	var totalDuration time.Duration
	for _, step := range l.steps {
		statusIcon := "✓"
		switch step.Status {
		case "error":
			statusIcon = "✗"
			errCount++
		case "warn":
			statusIcon = "△"
			warnCount++
		}
		fmt.Fprintf(l.logger.Writer(), "  %s %-20s  %-10s  %v\n",
			statusIcon, step.Name, step.Duration.Round(time.Millisecond), step.Status)
		totalDuration += step.Duration
	}
	fmt.Fprintf(l.logger.Writer(), "------------------------------\n")
	fmt.Fprintf(l.logger.Writer(), "  总计: %d 步骤, %d 错误, %d 警告, 耗时: %v\n\n",
		len(l.steps), errCount, warnCount, totalDuration.Round(time.Millisecond))
}

// HasErrors 检查是否有步骤失败
func (l *StepLogger) HasErrors() bool {
	for _, step := range l.steps {
		if step.Status == "error" {
			return true
		}
	}
	return false
}
