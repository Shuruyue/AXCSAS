# AXCSAS 專案總結分析報告
**(Project Analysis Report: Electroplated Copper Microstructure Evolution)**

本報告彙整了 AXCSAS 專案針對 **52 個電鍍銅樣品 (0~24小時自退火)** 的全方位分析結果。
分析涵蓋三大維度：織構 (Texture)、缺陷 (Defect)、與微觀結構 (Microstructure)。

---

## 1. 織構演變 (Texture Evolution)
**核心問題：添加劑如何改變銅的生長取向？**

### 圖表分析：[texture_TC_by_concentration.png](file:///d:/Shuru/Git%20Project/AXCSAS/outputs/plots/texture_analysis/texture_TC_by_concentration.png)
*   **觀察**：
    *   **0mL (無添加劑)**：(220) 訊號較強，呈現典型的生長模式。
    *   **18mL (強添加劑)**：(111) 訊號逐漸占主導地位 (TC > 1)。
*   **物理機制**：
    *   添加劑 (Leveler/SPS) 抑制了 (220) 和 (200) 的側向生長，迫使晶粒以最密排面 (111) 向上堆疊。
    *   這對於半導體銅製程至關重要，因為 (111) 取向具有最強的**抗電遷移 (Electromigration Resistance)** 能力。

---

## 2. 缺陷動力學 (Defect Kinetics)
**核心問題：強化的來源是什麼？它們如何隨時間消失？**

### A. 堆垛層錯：[stacking_fault_evolution.png](file:///d:/Shuru/Git%20Project/AXCSAS/outputs/plots/defect_analysis/stacking_fault_evolution.png)
*   **觀察**：
    *   **初始狀態 (0h)**：18mL 樣品的層錯機率高達 **0.4%** (每 50nm 一個缺陷)。這比 0mL 樣品高出許多。
    *   **演變 (0-10h)**：層錯機率呈現**S型下降**曲線，這是典型的 Johnson-Mehl-Avrami-Kolmogorov (JMAK) 再結晶動力學特徵。
    *   **終態 (24h)**：缺陷幾乎完全消失，回復到接近 0% 的完美晶體狀態。

### B. 殘餘應力：[stress_evolution.png](file:///d:/Shuru/Git%20Project/AXCSAS/outputs/plots/defect_analysis/stress_evolution.png)
*   **應力與層錯的解耦合 (Decoupling)**：
    *   **(111) 曲線 (~1000 MPa)**：數值異常偏高。經分析確認，這是由於 **堆垛層錯 (Stacking Faults)** 將衍射峰向高角度推移，疊加了虛假的拉伸應力分量。
    *   **(220) 曲線 (~250 MPa)**：**(重要)** 由於 (220) 平面不受層錯偏移影響，此曲線代表了樣品的**真實殘留應力**。
*   **演變趨勢**：
    *   應力隨時間呈現拉伸建立 (Tensile Buildup)，符合晶粒生長緻密化 (Densification) 模型。
    *   真實應力維持在 ~250 MPa，這是典型的添加劑強化電鍍銅特徵 (非超硬奈米孿晶，但具有高密度缺陷)。

---

## 3. 微觀應變 (Microstructure / Williamson-Hall)
**核心問題：晶格內部有多「亂」？**

### 圖表分析：[microstrain_evolution.png](file:///d:/Shuru/Git%20Project/AXCSAS/outputs/plots/wh_analysis/microstrain_evolution.png)
*   **觀察**：
    *   高濃度樣品 (18mL) 的微觀應變 (Microstrain) 初始值極高。
    *   隨著時間推移，應變急劇下降。
*   **結論**：
    *   這一結果與 Stack Fault 結果互為佐證。Stacking Fault 是微觀應變的主要貢獻者。
    *   應變的釋放對應著電阻率 (Resistivity) 的下降。

---

## 4. 綜合結論 (Overall Experiment Conclusion)

本研究成功利用非破壞性 XRD 技術，完整重建了電鍍銅的**自退火 (Self-annealing)** 過程：

1.  **添加劑效應**：成功誘導出高密度平面缺陷 ($\alpha \approx 0.4\%$) 和 (111) 優選取向。
2.  **物理模型修正**：釐清了異常高的 (111) 應力訊號實為層錯效應，並利用 (220) 訊號確立了真實應力水平 (~250 MPa)。這解答了樣品並非奈米孿晶但訊號異常的謎題。
3.  **動力學機制**：在室溫下，系統透過消除層錯與釋放應力來降低內能，導致了微觀結構的劇烈重組。

這套數據具有極高的內部一致性 (Self-consistency)，從織構、應力到層錯，所有指標都指向同一個物理故事。
