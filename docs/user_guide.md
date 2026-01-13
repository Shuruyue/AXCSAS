# AXCSAS 使用指南

## 快速開始

### 安裝

```bash
pip install -e .
```

### 基本使用

```bash
# 分析 XRD 資料
axcsas analyze data/raw/sample.txt -o outputs/

# 產生報告
axcsas report outputs/results.csv

# 查看幫助
axcsas --help
```

---

## 模組使用

### 1. Scherrer 晶粒尺寸分析

```python
from axcsas.methods import ScherrerCalculator

calc = ScherrerCalculator(wavelength=1.54056)
result = calc.calculate(
    two_theta=43.3,
    fwhm_observed=0.25,
    fwhm_instrumental=0.08,
    hkl=(1, 1, 1)
)

print(f"晶粒尺寸: {result.size_nm:.1f} nm")
print(f"有效性: {result.validity_flag}")
```

### 2. Williamson-Hall 分析

```python
from axcsas.methods import WilliamsonHallAnalyzer
import numpy as np

wh = WilliamsonHallAnalyzer()
two_theta = np.array([43.3, 50.4, 74.1, 89.9])
fwhm = np.array([0.25, 0.27, 0.32, 0.38])

result = wh.analyze(two_theta, fwhm)

print(f"晶粒尺寸: {result.crystallite_size_nm:.1f} nm")
print(f"微應變: {result.microstrain:.2e}")
print(f"R²: {result.r_squared:.3f}")
```

### 3. 織構分析

```python
from axcsas.methods import TextureAnalyzer

analyzer = TextureAnalyzer()
intensities = {
    (1, 1, 1): 100,
    (2, 0, 0): 35,
    (2, 2, 0): 15,
}

result = analyzer.analyze(intensities)

print(f"優選方向: {result.dominant_hkl}")
print(f"織構係數: {result.tc_values}")
```

### 4. 完整分析流程

```python
from axcsas.analysis import AXCSASPipeline

pipeline = AXCSASPipeline()
result = pipeline.analyze("data/raw/sample.txt")

print(result.report)
```

---

## 輸出結構

```
outputs/
├── analysis/       # CSV 分析結果
├── plots/          # 視覺化圖表
│   ├── fwhm/
│   ├── scherrer/
│   ├── texture/
│   └── williamson_hall/
├── reports/        # 文字報告
└── calibration/    # 校準資料
```

---

## 配置

編輯 `config.yaml` 調整分析參數：

```yaml
instrument:
  wavelength: 1.54056      # Cu Kα1
  caglioti_U: null         # 需校準
  caglioti_V: null
  caglioti_W: null

analysis:
  scherrer_k: 0.89
  wh_r2_threshold: 0.70

output:
  format: csv
  include_plots: true
```

---

## 常見問題

### Q: 儀器參數未校準怎麼辦？
A: 使用 LaB6 或 Si 標準樣品執行校準：
```bash
axcsas calibrate data/standards/LaB6.txt
```

### Q: 為什麼 Scherrer 結果標記為 UNRELIABLE？
A: 當計算的晶粒尺寸超出 5-200 nm 範圍時會標記。可能原因：
- FWHM 太小（接近儀器極限）
- FWHM 太大（奈米級以下晶粒）

### Q: W-H R² 很低怎麼辦？
A: 低 R² 通常表示：
- 數據點不足（需 ≥3 個峰）
- 存在各向異性應變
- 峰擬合不準確
