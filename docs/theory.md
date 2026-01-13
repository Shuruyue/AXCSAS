# AXCSAS 理論基礎

## 概述

AXCSAS (Advanced XRD Crystallite Size Analysis System) 基於 X 射線繞射 (XRD) 理論，實現銅電鍍薄膜的微結構分析。本文檔詳述各模組使用的數學公式與運算邏輯。

---

## 1. Scherrer 晶粒尺寸分析

### 1.1 基本公式

Scherrer 方程式用於從 XRD 峰寬估算晶粒尺寸：

$$
D = \frac{K \lambda}{\beta \cos\theta}
$$

**參數說明**：
| 符號 | 說明 | 單位 |
|------|------|------|
| $D$ | 晶粒尺寸 | nm |
| $K$ | Scherrer 常數 (形狀因子) | 無量綱 |
| $\lambda$ | X 射線波長 (Cu Kα₁ = 1.54056 Å) | Å |
| $\beta$ | 樣品線寬 (FWHM) | rad |
| $\theta$ | Bragg 角 (2θ/2) | rad |

### 1.2 Scherrer 常數 K

不同晶粒形狀對應不同 K 值：

| 形狀 | K 值 | 適用 hkl |
|------|------|----------|
| 球形 | 0.89 | 通用 |
| 立方體 (100) | 0.89 | (100), (200) |
| 立方體 (110) | 1.00 | (110), (220) |
| 立方體 (111) | 1.155 | (111) |

**實作邏輯**（`methods/scherrer.py`）：
```python
K_VALUES = {
    (1, 0, 0): 0.89,
    (1, 1, 0): 1.00,
    (1, 1, 1): 1.155,
    (2, 0, 0): 0.89,
    (2, 2, 0): 1.00,
    (3, 1, 1): 1.00,
}
```

### 1.3 儀器展寬校正

觀測線寬包含儀器展寬與樣品展寬：

$$
\beta_{sample}^2 = \beta_{obs}^2 - \beta_{inst}^2
$$

儀器展寬使用 **Caglioti 方程式** 計算：

$$
\beta_{inst}^2 = U \tan^2\theta + V \tan\theta + W
$$

**有效性判定**：
- $\beta_{obs} / \beta_{inst} \geq 1.2$：結果可靠
- $\beta_{obs} / \beta_{inst} < 1.2$：接近儀器極限，標記 WARNING
- $\beta_{obs} < \beta_{inst}$：無法計算，標記 ERROR

### 1.4 有效性旗標

| 旗標 | 條件 | 說明 |
|------|------|------|
| VALID | D ∈ [5, 200] nm | 結果可信 |
| WARNING | D ∈ [2, 5) 或 (200, 500] nm | 邊界值 |
| UNRELIABLE | D < 2 或 D > 500 nm | 超出適用範圍 |
| ERROR | 計算失敗 | 輸入無效 |

---

## 2. Williamson-Hall 微應變分析

### 2.1 基本原理

W-H 方法同時考慮晶粒尺寸與微應變對線寬的貢獻：

$$
\beta \cos\theta = \frac{K\lambda}{D} + 4\varepsilon \sin\theta
$$

線性化為 W-H 圖：
- **Y 軸**：$\beta \cos\theta$
- **X 軸**：$4 \sin\theta$
- **斜率**：微應變 $\varepsilon$
- **截距**：$K\lambda / D$

### 2.2 計算流程

1. **單位轉換**：
   $$\beta_{rad} = \beta_{deg} \times \frac{\pi}{180}$$

2. **座標計算**：
   $$x_i = 4 \sin(\theta_i), \quad y_i = \beta_i \cos(\theta_i)$$

3. **線性回歸**：
   $$y = mx + c$$
   
4. **參數提取**：
   $$\varepsilon = m, \quad D = \frac{K\lambda}{c}$$

### 2.3 品質評估

使用決定係數 $R^2$ 評估擬合品質：

| $R^2$ | 品質 | 說明 |
|-------|------|------|
| ≥ 0.90 | EXCELLENT | 線性關係良好 |
| 0.70 - 0.89 | GOOD | 可接受 |
| 0.50 - 0.69 | MARGINAL | 需謹慎解讀 |
| < 0.50 | POOR | 不可靠 |

**最少峰數要求**：≥ 3 個峰（2 個峰標記 WARNING）

### 2.4 各向異性診斷

$$
\sigma_{residual} = \sqrt{\frac{\sum(y_i - \hat{y}_i)^2}{n-2}}
$$

若殘差 > 閾值，可能存在各向異性應變。

---

## 3. 織構係數分析

### 3.1 Harris 織構係數

定義：

$$
TC_{hkl} = \frac{I_{hkl} / I_{hkl}^0}{\frac{1}{n} \sum_{i=1}^{n} I_i / I_i^0}
$$

**參數說明**：
| 符號 | 說明 |
|------|------|
| $TC_{hkl}$ | (hkl) 方向的織構係數 |
| $I_{hkl}$ | 實測 (hkl) 峰強度 |
| $I_{hkl}^0$ | JCPDS 標準強度 |
| $n$ | 參與計算的峰數 |

### 3.2 JCPDS 標準強度 (Cu, 04-0836)

| hkl | 2θ (°) | 標準強度 |
|-----|--------|----------|
| (111) | 43.30 | 100 |
| (200) | 50.43 | 46 |
| (220) | 74.13 | 20 |
| (311) | 89.93 | 17 |

### 3.3 守恆性質

$$
\sum_{i=1}^{n} TC_i = n
$$

### 3.4 織構程度

$$
\sigma = \sqrt{\frac{\sum(TC_i - 1)^2}{n}}
$$

| σ 值 | 判定 |
|------|------|
| < 0.1 | 隨機織構 |
| 0.1 - 0.5 | 輕微優選 |
| > 0.5 | 強烈優選 |

### 3.5 優選方向判定

$$
TC_{max} > 1 + 2\sigma \Rightarrow \text{存在優選方向}
$$

---

## 4. 缺陷與應力分析

### 4.1 堆垛層錯 (Stacking Fault)

**Warren 公式**：

$$
\alpha = \frac{\Delta 2\theta_{(111)-(200)}}{G}
$$

其中：
- $\Delta 2\theta = 2\theta_{exp} - 2\theta_{std}$
- $G = -0.5$ (Warren G 係數)
- 標準峰間距：$2\theta_{(200)} - 2\theta_{(111)} = 7.136°$

**嚴重程度分類**：
| 條件 | 分類 |
|------|------|
| α < 0.5% | NORMAL |
| 0.5% ≤ α < 2% | ELEVATED |
| 2% ≤ α < 5% | HIGH |
| α ≥ 5% | CRITICAL |

**SPS 預警**：峰間距 < 7.0° 觸發 Self-annealing Probability Signal

### 4.2 晶格常數計算

Bragg 定律：

$$
n\lambda = 2d \sin\theta
$$

面間距：

$$
d_{hkl} = \frac{a}{\sqrt{h^2 + k^2 + l^2}}
$$

晶格常數：

$$
a = d_{hkl} \times \sqrt{h^2 + k^2 + l^2}
$$

**Cu 標準晶格常數**：$a_0 = 3.6150$ Å

### 4.3 晶格狀態判定

$$
\delta = \frac{a - a_0}{a_0} \times 100\%
$$

| δ (%) | 狀態 |
|-------|------|
| < 0.05 | NORMAL |
| 0.05 - 0.20 | MINOR_EXPANSION |
| > 0.20 | EXPANSION |
| < -0.05 | CONTRACTION |

### 4.4 自退火狀態

銅電鍍膜的自退火時間表：

| 時間 | 狀態 |
|------|------|
| < 1 小時 | AS_DEPOSITED (鍍態) |
| 1-24 小時 | PARTIAL (部分退火) |
| 1-7 天 | ANNEALED (退火完成) |
| > 7 天 | STABLE (穩定) |

---

## 5. 峰擬合

### 5.1 Pseudo-Voigt 函數

$$
PV(x) = \eta \cdot L(x) + (1-\eta) \cdot G(x)
$$

**Lorentzian**：
$$
L(x) = \frac{I_0}{1 + 4\left(\frac{x-x_0}{\Gamma}\right)^2}
$$

**Gaussian**：
$$
G(x) = I_0 \exp\left[-4\ln 2 \left(\frac{x-x_0}{\Gamma}\right)^2\right]
$$

| 參數 | 說明 |
|------|------|
| $\eta$ | 混合參數 (0=純 Gaussian, 1=純 Lorentzian) |
| $I_0$ | 峰強度 |
| $x_0$ | 峰位置 |
| $\Gamma$ | FWHM |

### 5.2 背景擬合

線性背景：$B(x) = ax + b$

---

## 參考文獻

1. Langford, J. I., & Wilson, A. J. C. (1978). *Scherrer after sixty years*. J. Appl. Cryst., 11, 102-113.
2. Williamson, G. K., & Hall, W. H. (1953). *X-ray line broadening from filed aluminium and wolfram*. Acta Metallurgica, 1, 22-31.
3. Warren, B. E. (1969). *X-ray Diffraction*. Addison-Wesley.
4. JCPDS Card No. 04-0836 (Copper).
