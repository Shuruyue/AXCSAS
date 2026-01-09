# AXCSAS Theoretical Background

本文件說明 AXCSAS 系統中使用的物理模型與數學方程式。

---

## 1. X-ray Diffraction Line Broadening

XRD 繞射峰展寬來自兩個主要來源：

### 1.1 Instrumental Broadening (儀器展寬)

由光學元件、X-ray 源幾何形狀等儀器因素造成。使用 **Caglioti 方程式**描述：

$$FWHM_{inst}^2 = U \tan^2\theta + V \tan\theta + W$$

其中 U, V, W 由標準樣品（如 NIST SRM 660c LaB₆）校正取得。

### 1.2 Sample Broadening (樣品展寬)

由晶粒尺寸與微觀應變造成：

- **Size broadening**：晶粒尺寸越小，峰越寬
- **Strain broadening**：晶格應變導致 d-spacing 變化

---

## 2. Scherrer Equation

計算晶粒尺寸的基本公式：

$$D = \frac{K\lambda}{\beta\cos\theta}$$

| 參數 | 說明 | 單位 |
|------|------|------|
| D | 晶粒尺寸 | Å 或 nm |
| K | 形狀因子 | 無因次 |
| λ | X-ray 波長 | Å |
| β | 積分峰寬（弧度） | rad |
| θ | Bragg 角 | rad |

### 形狀因子 K

| 晶粒形狀 | K 值 | 來源 |
|----------|------|------|
| 球形 | 0.89 | Langford & Wilson (1978) |
| 立方體 | 0.94 | Langford & Wilson (1978) |

---

## 3. Williamson-Hall Analysis

同時分離尺寸與應變貢獻：

$$\beta\cos\theta = \frac{K\lambda}{D} + 4\varepsilon\sin\theta$$

將 $\beta\cos\theta$ 對 $\sin\theta$ 作圖：
- **截距** = $K\lambda/D$ → 計算晶粒尺寸
- **斜率** = $4\varepsilon$ → 計算微觀應變

---

## 4. Pseudo-Voigt Profile

真實 XRD 峰型為 Gaussian 與 Lorentzian 的組合：

$$I(2\theta) = I_0 \cdot [\eta \cdot L(2\theta) + (1-\eta) \cdot G(2\theta)]$$

| 分量 | 物理來源 |
|------|----------|
| Lorentzian (L) | 晶粒尺寸展寬 |
| Gaussian (G) | 微觀應變 + 儀器效應 |
| η | 混合參數 (0-1) |

---

## 5. Harris Texture Coefficient

量化優選取向程度：

$$TC(hkl) = \frac{I(hkl)/I_0(hkl)}{\frac{1}{n}\sum I(hkl)/I_0(hkl)}$$

**解釋：**
- TC = 1：隨機取向（粉末平均）
- TC > 1：此方向為優選取向
- TC < 1：此方向被抑制

---

## 6. Instrumental Broadening Correction

對於 Pseudo-Voigt 峰型，使用幾何近似法：

$$\beta_{sample} = \beta_{obs} - \frac{\beta_{inst}^2}{\beta_{obs}}$$

**有效性條件：**
- $\beta_{obs} > 1.2 \times \beta_{inst}$
- 在此範圍內誤差 < 1%

---

## References

1. Langford, J. I., & Wilson, A. J. C. (1978). *Scherrer after sixty years*. J. Appl. Cryst., 11, 102-113.
2. Williamson, G. K., & Hall, W. H. (1953). *X-ray line broadening*. Acta Metallurgica, 1(1), 22-31.
3. Caglioti, G., Paoletti, A., & Ricci, F. P. (1958). *Choice of collimators for a crystal spectrometer*. Nuclear Instruments, 3(4), 223-228.
