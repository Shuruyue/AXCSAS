# AXCSAS 
## Advanced XRD Crystallite Size Analysis System

基於偽沃伊特卷積與儀器校正的自動化晶粒尺寸計算系統

---

## 📋 專案簡介

AXCSAS 是一套自動化演算法系統，解決傳統手動計算 Scherrer Size 時面臨的：
- 基線選取主觀誤差
- 儀器展寬扣除不精確
- 峰型擬合函數選擇錯誤

### 核心特色
- **Pseudo-Voigt 全譜擬合**：採用學術界公認最接近真實晶體繞射行為的數學描述
- **Caglioti 儀器校正**：利用 NIST SRM 660c (LaB₆) 進行全角度儀器寬度校正
- **高精度預測**：適用於 2-100 nm 範圍內晶粒尺寸

---

## 🚀 快速開始

### 安裝依賴
```bash
pip install -r requirements.txt
```

### 儀器校正
```bash
python scripts/calibrate_instrument.py --standard data/standards/LaB6_SRM660c.xy
```

### 樣品分析
```bash
python scripts/analyze_sample.py --input data/raw/sample.xy --output outputs/results/
```

### 批次分析
```bash
python scripts/batch_analysis.py --input-dir data/raw/202511/ --output-dir outputs/results/
```

---

## 📁 專案結構

```
AXCSAS/
├── config.yaml              # 全域設定檔
├── data/                    # 數據目錄
│   ├── raw/                 # 原始 XRD 數據
│   ├── standards/           # NIST 標準品數據
│   └── processed/           # 預處理後數據
├── src/                     # 核心程式碼
│   ├── preprocessing/       # 數據預處理模組
│   ├── fitting/             # 峰值擬合核心
│   ├── physics/             # 物理計算核心
│   ├── validation/          # 誤差分析與驗證
│   └── utils/               # 工具函式
├── scripts/                 # 執行腳本
├── outputs/                 # 輸出目錄
├── tests/                   # 單元測試
└── docs/                    # 文件
```

---

## 📐 理論基礎

### Pseudo-Voigt 峰型函數
$$I(2\theta) = I_0 \cdot [ \eta L(2\theta) + (1-\eta) G(2\theta) ] + Background$$

### Caglioti 方程式
$$FWHM_{inst}^2 = U \tan^2\theta + V \tan\theta + W$$

### Scherrer 方程式
$$D = \frac{K \lambda}{\beta \cos\theta}$$

---

## 📚 參考文獻

- Langford, J. I., & Wilson, A. J. C. (1978). *Scherrer after sixty years*. J. Appl. Cryst., 11, 102-113.
- NIST Standard Reference Material 660c (LaB₆)

---

## 📄 授權

MIT License
