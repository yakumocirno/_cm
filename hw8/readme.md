# [AI對話](https://chatgpt.com/share/69597be2-a23c-8005-84a3-6e6ed20a336e)
好，下面這一份是**「可直接複製貼上、不含任何額外標記」的 `README.md` 原文**。
你可以 **整段全選 → Ctrl+C → 貼到 README.md**，可直接用於作業或專案。

---

# 資訊理論與編碼

## 機率、熵、交叉熵、KL 散度與 (7,4) Hamming Code 實作說明

本專案以 **Python** 實作資訊理論與通道編碼的核心概念，並透過程式驗證其數學性質。內容包含極小機率的數值計算、對數穩定化、熵與交叉熵、KL 散度、互資訊，以及 **(7,4) 漢明碼** 的編碼與單錯誤更正解碼。

---

## 1. 公平銅板連續 10000 次全正面的機率

### 數學原理

對一枚公平銅板
P(正面) = p = 0.5

連續丟擲 n = 10000 次全部為正面的機率為：

p^10000 = 2^(-10000)

這是一個極端微小的數值，遠小於 IEEE 754 double precision 能表示的最小正數，因此在實際計算中會發生 **浮點下溢（underflow）**。

### 程式設計

```python
def prob_all_heads(p: float, n: int) -> float:
    return p ** n
```

程式中同時提供 Monte Carlo 模擬版本，用來說明在有限次模擬下，幾乎不可能觀測到此事件。

---

## 2. 使用對數避免數值下溢

### 數學原理

利用對數性質：

log(p^n) = n log(p)

即使 p^n 無法直接表示，log(p^n) 仍可精確計算，這是資訊理論與機率模型中常用的技巧。

### 程式設計

```python
def log_prob_power(p: float, n: int) -> float:
    return n * math.log(p)
```

並示範轉換為：

* 自然對數
* log10
* log2

同時說明將對數轉回機率仍會發生下溢。

---

## 3. 資訊理論量的實作

### (1) 熵 Entropy

定義為：
H(p) = −∑ p_i log p_i

代表隨機變數的不確定性下界。

```python
def entropy(p: List[float], base=2):
```

---

### (2) 交叉熵 Cross Entropy

定義為：
H(p, q) = −∑ p_i log q_i

表示「使用分佈 q 來編碼來自分佈 p 的資料時，平均需要的碼長」。

```python
def cross_entropy(p, q, base=2):
```

---

### (3) KL 散度 Kullback–Leibler Divergence

定義為：
D_KL(p || q) = ∑ p_i log(p_i / q_i)

衡量用 q 近似 p 所造成的資訊損失，且永遠 ≥ 0。

```python
def kl_divergence(p, q, base=2):
```

---

### (4) 互資訊 Mutual Information

定義為：
I(X;Y) = ∑ p(x,y) log( p(x,y) / (p(x)p(y)) )

表示兩個隨機變數之間共享的資訊量。

```python
def mutual_information(joint, base=2):
```

---

## 4. 驗證交叉熵不等式

### 理論關係

交叉熵可分解為：
H(p, q) = H(p) + D_KL(p || q)

因為 KL 散度永遠 ≥ 0，因此：
H(p, q) ≥ H(p, p)

且當 q ≠ p 時嚴格大於。

### 程式驗證

```python
verify_cross_entropy_inequality(p_dist, q_dist)
```

輸出數值結果，驗證理論不等式方向正確。

---

## 5. (7,4) Hamming Code 編碼與解碼

### 編碼背景

(7,4) 漢明碼使用：

* 4 個資料位元
* 3 個 parity 位元
  可更正 **單一位元錯誤（1-bit error correction）**。

### 位元配置（1-based index）

位置：1  2  3  4  5  6  7
類型：p1 p2 d1 p4 d2 d3 d4

### EVEN parity 設計

p1 = d1 ⊕ d2 ⊕ d4
p2 = d1 ⊕ d3 ⊕ d4
p4 = d2 ⊕ d3 ⊕ d4

### 解碼（Syndrome）

s1 = XOR(1,3,5,7)
s2 = XOR(2,3,6,7)
s4 = XOR(4,5,6,7)

錯誤位置 = s1 + 2·s2 + 4·s4

### 程式實作

```python
hamming74_encode(data4)
hamming74_decode(code7)
```

示範流程包含：

* 編碼
* 人為製造 1-bit error
* 成功定位並修正錯誤

---

## 6. 香農相關定理說明

### (1) 香農通道編碼定理

只要資料傳輸速率 R 小於通道容量 C，就存在一種編碼方式，使得錯誤率可以趨近於 0。

此定理說明「**可靠通訊在理論上是可行的**」。

---

### (2) 香農–哈特利定理

對頻寬為 B、訊雜比為 SNR 的通道，其容量為：

C = B · log2(1 + SNR)

此公式給出實體通道可達的最大資訊傳輸速率。

---

## 總結

本專案完成以下目標：

* 展示極小機率的數值計算與對數技巧
* 實作熵、交叉熵、KL 散度與互資訊
* 驗證交叉熵與 KL 散度的理論不等式
* 實作 (7,4) 漢明碼並成功進行單錯誤更正
* 說明香農定理在通訊理論中的核心意義

本程式完整串聯 **數學理論 → 演算法 → 程式實作 → 數值驗證**。

---

如果你需要
✔ 精簡成「作業版」
✔ 改寫成「數學報告格式（含公式編號）」
✔ 或補上「線性代數視角的漢明碼推導」

直接跟我說即可。
