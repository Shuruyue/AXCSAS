# AXCSAS User Guide

AXCSAS (Advanced XRD Crystallite Size Analysis System) 使用指南。

---

## Quick Start

### 1. 安裝依賴

```bash
cd AXCSAS
pip install -r requirements.txt
```

### 2. 準備數據

將 XRD 數據檔案放入 `data/raw/` 目錄。

支援格式：`.xy`, `.csv`, `.txt` (Bruker)

### 3. 分析樣品

```bash
# 單一樣品
python scripts/analyze_sample.py -i data/raw/sample.xy -o outputs/results/

# 批次分析
python scripts/batch_analysis.py -i data/raw/ -o outputs/results/
```

### 4. 產生報告

```bash
python scripts/generate_report.py -i outputs/results/ -o outputs/reports/
```

---

## Command Reference

### analyze_sample.py

```
用法: python scripts/analyze_sample.py [OPTIONS]

選項:
  -i, --input     輸入 XRD 檔案路徑 (必填)
  -o, --output    輸出目錄 (預設: outputs/results)
  -c, --config    設定檔路徑 (預設: config.yaml)
  --no-plots      不產生圖表
```

### batch_analysis.py

```
用法: python scripts/batch_analysis.py [OPTIONS]

選項:
  -i, --input-dir   輸入目錄 (必填)
  -o, --output-dir  輸出目錄 (預設: outputs/results)
  -c, --config      設定檔路徑 (預設: config.yaml)
```

### generate_report.py

```
用法: python scripts/generate_report.py [OPTIONS]

選項:
  -i, --input   分析結果目錄 (必填)
  -o, --output  報告輸出目錄 (預設: outputs/reports)
  -t, --title   報告標題
```

---

## Configuration

編輯 `config.yaml` 調整分析參數：

```yaml
# 預處理設定
preprocessing:
  smoothing:
    window_size: 11    # Savitzky-Golay 窗口 (奇數)
    poly_order: 3      # 多項式階數
  background:
    method: "chebyshev"
    poly_degree: 5

# 峰值擬合
fitting:
  function: "pseudo_voigt"
  peak_detection:
    min_height: 100    # 最小峰高
    min_distance: 0.5  # 最小峰間距 (度)

# 驗證閾值
validation:
  min_broadening_ratio: 1.2
  max_rwp: 10.0
```

---

## Output Files

分析後產生的檔案：

| 檔案 | 說明 |
|------|------|
| `*_peaks.csv` | 峰值擬合結果 |
| `*_summary.yaml` | 分析摘要 |
| `*_fit.png` | 擬合圖 |
| `batch_summary.csv` | 批次匯總 |
| `analysis_report.md` | 完整報告 |

---

## Troubleshooting

### 問題：Size > 200 nm 警告

**原因**：峰寬太小，接近儀器極限  
**解決**：確認數據品質，考慮儀器校正

### 問題：W-H R² 過低

**原因**：峰數不足或非均勻應變  
**解決**：確保至少 3 個峰，檢查峰位判定

### 問題：Texture TC 異常

**原因**：強度數據不正確  
**解決**：確認峰高讀取正確，檢查 hkl 指定

---

## API Reference

```python
from src.physics import ScherrerCalculator, WilliamsonHallAnalyzer
from src.validation import ErrorAnalyzer

# Scherrer 計算
calc = ScherrerCalculator(wavelength=1.54056, k_factor=0.89)
result = calc.calculate(two_theta=43.3, fwhm=0.25)
print(f"Size: {result.size_nm:.1f} nm")

# W-H 分析
wh = WilliamsonHallAnalyzer()
result = wh.analyze(two_theta_array, fwhm_array)
print(f"Size: {result.crystallite_size_nm:.1f} nm")
print(f"Strain: {result.microstrain:.2e}")
```

---

## Support

問題回報：建立 Issue 或聯繫開發團隊
