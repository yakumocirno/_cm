# -*- coding: utf-8 -*-
"""
需求清單：
1) 公平銅板 p=0.5，連續投擲 10000 次全正面的機率：p^10000
2) 用 log(p^n) = n log(p) 計算 log(0.5^10000)
3) 計算：熵 H、交叉熵 CE、KL 散度、互資訊 I(X;Y)
4) 驗證：cross_entropy(p,p) < cross_entropy(p,q) 當 q != p 時（注意方向！）
5) (7,4) Hamming code 編碼與解碼範式（單錯誤更正）
6) 說明：夏農通道編碼定理、夏農-哈特利定理
"""

import math
import random
from typing import List, Tuple

# =========================
# 1) p^10000 與模擬
# =========================

def prob_all_heads(p: float, n: int) -> float:
    """直接算 p^n（可能非常小到變成 0.0，屬正常浮點下溢）"""
    return p ** n

def simulate_all_heads(p: float, n: int, trials: int = 1000) -> float:
    """
    蒙地卡羅估計：在 trials 次實驗中，丟 n 次都正面的比例。
    對 n=10000 幾乎永遠都是 0（因為真實機率 2^-10000 太小）。
    """
    cnt = 0
    for _ in range(trials):
        ok = True
        for _ in range(n):
            if random.random() >= p:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt / trials

# =========================
# 2) log(p^n) = n log(p)
# =========================

def log_prob_power(p: float, n: int) -> float:
    """回傳 log(p^n) ，用 n*log(p)（數值更穩定）"""
    if p <= 0:
        return float("-inf")
    return n * math.log(p)

# =========================
# 3) 資訊理論：熵、交叉熵、KL、互資訊
# =========================

EPS = 1e-15

def _normalize(dist: List[float]) -> List[float]:
    s = sum(dist)
    if s <= 0:
        raise ValueError("Distribution sum must be positive.")
    return [x / s for x in dist]

def entropy(p: List[float], base: float = 2.0) -> float:
    """H(p) = - sum p_i log p_i"""
    p = _normalize(p)
    logb = math.log(base)
    h = 0.0
    for pi in p:
        if pi > 0:
            h -= pi * (math.log(pi) / logb)
    return h

def cross_entropy(p: List[float], q: List[float], base: float = 2.0) -> float:
    """H(p,q) = - sum p_i log q_i"""
    p = _normalize(p)
    q = _normalize(q)
    logb = math.log(base)
    ce = 0.0
    for pi, qi in zip(p, q):
        # 若 q_i=0 且 p_i>0 交叉熵為 +inf；用極小值避免崩潰但仍顯示很大
        qi_safe = max(qi, EPS)
        if pi > 0:
            ce -= pi * (math.log(qi_safe) / logb)
    return ce

def kl_divergence(p: List[float], q: List[float], base: float = 2.0) -> float:
    """D_KL(p||q) = sum p_i log(p_i/q_i)"""
    p = _normalize(p)
    q = _normalize(q)
    logb = math.log(base)
    kl = 0.0
    for pi, qi in zip(p, q):
        if pi > 0:
            qi_safe = max(qi, EPS)
            kl += pi * (math.log(pi / qi_safe) / logb)
    return kl

def mutual_information(joint: List[List[float]], base: float = 2.0) -> float:
    """
    I(X;Y) = sum_{x,y} p(x,y) log( p(x,y) / (p(x)p(y)) )
    joint: 二維聯合分佈矩陣
    """
    # normalize joint
    total = sum(sum(row) for row in joint)
    if total <= 0:
        raise ValueError("Joint distribution must sum to positive.")
    Pxy = [[v / total for v in row] for row in joint]
    px = [sum(row) for row in Pxy]
    py = [sum(Pxy[i][j] for i in range(len(Pxy))) for j in range(len(Pxy[0]))]

    logb = math.log(base)
    I = 0.0
    for i in range(len(Pxy)):
        for j in range(len(Pxy[0])):
            pxy = Pxy[i][j]
            if pxy > 0:
                denom = max(px[i] * py[j], EPS)
                I += pxy * (math.log(pxy / denom) / logb)
    return I

def verify_cross_entropy_inequality(p: List[float], q: List[float]) -> None:
    """
    正確的不等式是：
      H(p, q) >= H(p, p)  (且 q!=p 時嚴格大於；在 q 有零且 p>0 會到 +inf)
    因為 H(p,q) = H(p) + KL(p||q)
          H(p,p) = H(p) + KL(p||p)=H(p)
          KL >= 0
    """
    ce_pp = cross_entropy(p, p)
    ce_pq = cross_entropy(p, q)
    print("cross_entropy(p,p) =", ce_pp)
    print("cross_entropy(p,q) =", ce_pq)
    if abs(ce_pq - ce_pp) < 1e-12:
        print("結果：幾乎相等（可能 q≈p）")
    elif ce_pq > ce_pp:
        print("驗證成功：cross_entropy(p,q) > cross_entropy(p,p) 當 q != p")
    else:
        print("驗證失敗：你給的 q 可能有特殊情況或數值誤差太大")

# =========================
# 4) (7,4) Hamming code：編碼與解碼
# =========================
"""
(7,4) 漢明碼：用 4 個資料位元 d1 d2 d3 d4
放在位置 3,5,6,7
parity 位元 p1,p2,p4 放在位置 1,2,4（位置用 1-based）
位元排列： [1]p1 [2]p2 [3]d1 [4]p4 [5]d2 [6]d3 [7]d4

奇偶校驗（採 EVEN parity，讓每組 XOR=0）：
p1 覆蓋位置 1,3,5,7 -> p1 = d1 XOR d2 XOR d4
p2 覆蓋位置 2,3,6,7 -> p2 = d1 XOR d3 XOR d4
p4 覆蓋位置 4,5,6,7 -> p4 = d2 XOR d3 XOR d4

解碼：算 syndrome s1,s2,s4：
s1 = XOR(1,3,5,7)
s2 = XOR(2,3,6,7)
s4 = XOR(4,5,6,7)
錯誤位置 = s1 + 2*s2 + 4*s4（這是二進位位置）
若不為 0，就翻轉該位置修正（可更正 1-bit error）
"""

def _xor_bits(bits: List[int]) -> int:
    x = 0
    for b in bits:
        x ^= (b & 1)
    return x

def hamming74_encode(data4: List[int]) -> List[int]:
    """data4: [d1,d2,d3,d4] -> code7: [c1..c7]"""
    if len(data4) != 4:
        raise ValueError("data4 must have length 4.")
    d1, d2, d3, d4 = [b & 1 for b in data4]
    p1 = d1 ^ d2 ^ d4
    p2 = d1 ^ d3 ^ d4
    p4 = d2 ^ d3 ^ d4
    code = [p1, p2, d1, p4, d2, d3, d4]
    return code

def hamming74_decode(code7: List[int]) -> Tuple[List[int], int, List[int]]:
    """
    回傳 (data4, error_pos, corrected_code7)
    error_pos: 0 表示無錯；1..7 表示更正了那個位置
    """
    if len(code7) != 7:
        raise ValueError("code7 must have length 7.")
    c = [b & 1 for b in code7]

    # syndrome bits (EVEN parity -> XOR should be 0)
    s1 = _xor_bits([c[0], c[2], c[4], c[6]])  # 1,3,5,7
    s2 = _xor_bits([c[1], c[2], c[5], c[6]])  # 2,3,6,7
    s4 = _xor_bits([c[3], c[4], c[5], c[6]])  # 4,5,6,7

    error_pos = s1 + 2 * s2 + 4 * s4  # 1..7
    corrected = c[:]
    if error_pos != 0:
        idx = error_pos - 1
        corrected[idx] ^= 1  # flip

    # extract data bits from corrected positions 3,5,6,7
    d1 = corrected[2]
    d2 = corrected[4]
    d3 = corrected[5]
    d4 = corrected[6]
    return [d1, d2, d3, d4], error_pos, corrected

# =========================
# Demo / 主程式
# =========================

def main():
    print("=== 1) 公平銅板連續 10000 次全正面機率 ===")
    p = 0.5
    n = 10000
    direct = prob_all_heads(p, n)
    print("p^10000 =", direct)
    # 這通常會直接下溢變成 0.0（因為太小）
    print("（浮點可能顯示 0.0 是正常的數值下溢）")

    print("\n=== 2) 用 log(p^n)=n log(p) 計算 log(0.5^10000) ===")
    lp = log_prob_power(p, n)
    print("log(0.5^10000) (natural log) =", lp)
    print("轉成 log10 =", lp / math.log(10))
    print("轉成 log2  =", lp / math.log(2))
    # 若你想要把它轉回機率（會下溢）
    print("exp(log) =", math.exp(lp))

    print("\n=== 3) 熵 / 交叉熵 / KL / 互資訊 範例 ===")
    p_dist = [0.1, 0.2, 0.7]
    q_dist = [0.2, 0.2, 0.6]
    print("H(p) =", entropy(p_dist))
    print("H(p,q) =", cross_entropy(p_dist, q_dist))
    print("KL(p||q) =", kl_divergence(p_dist, q_dist))
    # mutual information example joint
    joint = [
        [0.1, 0.2],
        [0.3, 0.4],
    ]
    print("I(X;Y) =", mutual_information(joint))

    print("\n=== 4) 驗證 cross_entropy(p,q) > cross_entropy(p,p) 當 q != p ===")
    verify_cross_entropy_inequality(p_dist, q_dist)

    print("\n=== 5) (7,4) Hamming code 編碼/解碼 demo ===")
    data = [1, 0, 1, 1]  # d1 d2 d3 d4
    code = hamming74_encode(data)
    print("data =", data)
    print("encoded code7 =", code)

    # introduce 1-bit error
    code_err = code[:]
    flip_pos = 6  # 1..7
    code_err[flip_pos - 1] ^= 1
    print("with 1-bit error at pos", flip_pos, "->", code_err)

    data_hat, errpos, corrected = hamming74_decode(code_err)
    print("decoder syndrome says error_pos =", errpos)
    print("corrected code7 =", corrected)
    print("decoded data4 =", data_hat)

    print("\n=== 6) 香農定理說明（看下方文字） ===")

if __name__ == "__main__":
    main()
