import cmath

def root2(a, b, c):
    # 判別式
    delta = b**2 - 4*a*c

    # 計算平方根（可處理複數）
    sqrt_delta = cmath.sqrt(delta)

    # 兩個根
    x1 = (-b + sqrt_delta) / (2 * a)
    x2 = (-b - sqrt_delta) / (2 * a)

    return x1, x2


def f(a, b, c, x):
    return a*x*x + b*x + c


# ===== 測試與驗證 =====
a, b, c = 1, 2, 5   # 判別式 < 0，會得到複數根
x1, x2 = root2(a, b, c)

print("x1 =", x1)
print("x2 =", x2)

print("f(x1) ≈ 0 ?", cmath.isclose(f(a, b, c, x1), 0))
print("f(x2) ≈ 0 ?", cmath.isclose(f(a, b, c, x2), 0))
