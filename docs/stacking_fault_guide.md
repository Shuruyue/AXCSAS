# 堆垛層錯 (Stacking Fault) 完全解析
**(From Fundamentals to Engineering Applications)**

您問到了最核心的關鍵。這不僅是學術理論，更是**半導體電鍍工程 (Damascene Process)** 中控制銅線性能的秘密武器。

---

## 21. 真相大白：為什麼 (111) 應力這麼高？(The "Fake" Stress)

您說：**「我的不是奈米孿晶，所以我才懷疑算錯。」**
這句話點醒了我。我重新檢視了您的數據，發現了一個**巨大的物理耦合現象**。

**結論先說：您的直覺是對的。那個 1000 MPa 是「虛」的。**
**真實的應力應該只有 200~250 MPa 左右 (看 (200) 和 (220) 峰)。**

### A. 兇手是誰？(Coupling Effect)
您的算法沒有錯，但是 **XRD 的物理特性** 開了一個玩笑：
1.  **拉伸應力 (Tensile Stress)** 會讓所有峰往高角度跑 (d 變小)。
2.  **堆垛層錯 (Stacking Faults)** 會有特殊偏移動作：
    *   **(111) 峰**：被用力**推向高角度** (看起來像超高拉伸應力)。
    *   **(200) 峰**：被用力**推向低角度** (抵消了拉伸應力)。

### B. 證據確鑿 (Evidence in Your Plot)
請回頭看您的應力演變圖：
*   **(111) 曲線**：高高掛在 **1000 MPa** $\rightarrow$ 這是「真實應力 + 層錯偏移」的疊加結果 (Double Count)。
*   **(200) 曲線**：只有 **200 MPa** $\rightarrow$ 這是「真實應力 - 層錯偏移」的抵消結果。
*   **(220) 曲線**：約 **250 MPa** $\rightarrow$ (220) 不受層錯偏移影響，它是**最誠實的證人**。

### C. 您該如何解釋？(The Correct Narrative)
所以，您的樣品**並不是** 1000 MPa 的超硬奈米孿晶銅 (您是對的)。
它是一個 **應力約 200 MPa 的普通電鍍銅**，但是含有 **0.4% 的高密度層錯**。

在論文中，您必須這樣澄清，這會顯示您非常懂 XRD：

> "The apparent high tensile stress observed in the (111) direction (~1000 MPa) is attributed to the **peak shift contribution from stacking faults**, which shifts the (111) peak to higher angles. The **(220) peak**, being unaffected by stacking faults, provides a more accurate estimation of the true residual stress, which is approximately **250 MPa**. This confirms that the distinct peak shifts are primarily defect-driven rather than purely elastic strain."
>
> (翻譯：在 (111) 方向觀察到的表面高拉伸應力 (~1000 MPa) 歸因於**堆垛層錯的峰位偏移貢獻**，其將 (111) 峰推向高角度。**(220) 峰**由於不受層錯影響，提供了更準確的真實殘留應力估算，約為 **250 MPa**。這證實了顯著的峰位偏移主要由缺陷驅動，而非純粹的彈性應變。)

**最終結論：**
*   算法沒錯 (公式是對的)。
*   數據沒錯 (峰位是真的)。
*   **是 (111) 被層錯「汙染」了。** 請以 (220) 的 250 MPa 作為您的真實應力數據。這樣就完全符合您「非奈米孿晶」的認知了。

---

## 20. 文獻鐵證 (Literature Verification)
*(如前述，略)*

## 19. 為什麼應力不會歸零？ (Why Not Zero?)
*(如前述，略)*

## 18. 目前的數據足以說明了嗎？ (Sufficiency & Gap Analysis)
*(如前述，略)*

## 17. 所以那些圖沒用？ (So, are they useless?)
*(如前述，略)*

## 16. 儀器限制與圖表缺陷 (Limitations & "Just XRD")
*(如前述，略)*

## 15. 現在可信度高嗎？ (Is the Credibility High Now?)
*(如前述，略)*

## 14. 楊氏模數驗證 (E=191 GPa? Is it correct?)
*(如前述，略)*

## 13. 這個 % 高嗎？ (Is 0.4% High?)
*(如前述，略)*

## 12. 兩個技術問題的解答 (Technical Q&A)
*(如前述，略)*

## 11. 萬一其他也有扭曲 (Stress) 呢？ (The Ultimate Defense)
*(如前述，略)*

## 10. 針對您疑慮的終極解答 (Final Answers to Your Doubts)
*(如前述，略)*
