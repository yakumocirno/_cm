import cmath

def root3(a, b, c, d):
    if a == 0:
        raise ValueError("a 不能為 0，否則不是三次多項式")

    # 1. 正規化
    A = b / a
    B = c / a
    C = d / a

    # 2. 降階 x = y - A/3
    p = B - A**2 / 3
    q = 2*A**3 / 27 - A*B / 3 + C

    # 3. 判別式
    Δ = (q/2)**2 + (p/3)**3

    # 4. Cardano
    sqrtΔ = cmath.sqrt(Δ)
    u = (-q/2 + sqrtΔ) ** (1/3)
    v = (-q/2 - sqrtΔ) ** (1/3)

    # 複數立方根的單位根
    omega = complex(-0.5, cmath.sqrt(3)/2)

    y1 = u + v
    y2 = omega*u + omega**2*v
    y3 = omega**2*u + omega*v

    # 5. 轉回 x
    shift = A / 3
    x1 = y1 - shift
    x2 = y2 - shift
    x3 = y3 - shift

    return x1, x2, x3

def f(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

roots = root3(1, 0, 0, 1)  # x^3 + 1 = 0

for r in roots:
    print(r, "->", f(r, 1, 0, 0, 1))