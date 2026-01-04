import cmath
import math

# =========================
# 參數設定
# =========================

L = 10.0      # x 積分範圍 [-L, L]
W = 10.0      # omega 積分範圍 [-W, W]
N = 2000      # x 取樣點數
M = 2000      # omega 取樣點數

dx = 2 * L / N
dw = 2 * W / M

# =========================
# 測試函數 f(x)
# =========================

def f(x):
    # 高斯函數（Fourier 轉換性質很好）
    return math.exp(-x * x)

# =========================
# 正傅立葉轉換（數值積分）
# F(ω) = ∫ f(x) e^{-i ω x} dx
# =========================

def dft(f):
    omegas = [(-W + j * dw) for j in range(M)]
    F = []

    for w in omegas:
        s = 0j
        for i in range(N):
            x = -L + i * dx
            s += f(x) * cmath.exp(-1j * w * x) * dx
        F.append(s)

    return omegas, F

# =========================
# 逆傅立葉轉換（數值積分）
# f(x) = (1/2π) ∫ F(ω) e^{i ω x} dω
# =========================

def idft(omegas, F):
    xs = [(-L + i * dx) for i in range(N)]
    f_rec = []

    for x in xs:
        s = 0j
        for w, Fw in zip(omegas, F):
            s += Fw * cmath.exp(1j * w * x) * dw
        f_rec.append((1 / (2 * math.pi)) * s)

    return xs, f_rec

# =========================
# 驗證：f -> F -> f
# =========================

omegas, F = dft(f)
xs, f_reconstructed = idft(omegas, F)

# =========================
# 驗證誤差
# =========================

max_error = 0.0
for x, fr in zip(xs, f_reconstructed):
    err = abs(fr.real - f(x))
    max_error = max(max_error, err)

print("最大誤差 =", max_error)
