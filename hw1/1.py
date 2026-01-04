# -*- coding: utf-8 -*-
import math

# -------------------------
# 數值積分：梯形法
# -------------------------
def integrate_trapezoid(f, a, b, n=4000):
    """
    近似積分 ∫_a^b f(t) dt (梯形法)
    n: 分割數，越大越準 (但太大會慢)
    """
    if a == b:
        return 0.0
    if b < a:
        return -integrate_trapezoid(f, b, a, n)

    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        s += f(a + i * h)
    return s * h

# -------------------------
# F(x) = ∫_0^x f(t) dt
# -------------------------
def F_of_x(f, x, n=4000):
    return integrate_trapezoid(f, 0.0, x, n)

# -------------------------
# 數值微分：中心差分
# -------------------------
def derivative_center(g, x, h=1e-4):
    return (g(x + h) - g(x - h)) / (2.0 * h)

# -------------------------
# 驗證：比較 d/dx ∫_0^x f(t)dt 與 f(x)
# -------------------------
def verify_ftc(f, xs, n_int=4000, h=1e-4):
    # 把 F 固定成一個函數（內部用數值積分）
    def G(x):
        return F_of_x(f, x, n=n_int)

    print("x\t\tF'(x)≈d/dx∫_0^x f(t)dt\tf(x)\t\tabs error")
    print("-" * 80)
    for x in xs:
        fp = derivative_center(G, x, h=h)
        fx = f(x)
        err = abs(fp - fx)
        print(f"{x: .5f}\t{fp: .10f}\t\t\t{fx: .10f}\t{err: .3e}")

# -------------------------
# 測試函數（你也可以自己換）
# -------------------------
def f1(t):
    return math.sin(t) + t*t

def f2(t):
    return math.exp(-t*t)

if __name__ == "__main__":
    xs = [-1.5, -1.0, -0.3, 0.0, 0.2, 0.7, 1.3, 2.0]

    print("\n=== 驗證 f(t) = sin(t) + t^2 ===")
    verify_ftc(f1, xs, n_int=6000, h=1e-4)

    print("\n=== 驗證 f(t) = exp(-t^2) ===")
    verify_ftc(f2, xs, n_int=6000, h=1e-4)
