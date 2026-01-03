def df(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2 * h)
def integral(f, a, b, n=10000):
    h = (b - a) / n
    s = 0.0
    for i in range(n):
        x1 = a + i * h
        x2 = x1 + h
        s += (f(x1) + f(x2)) * h / 2
    return s
def theorem1(f, x, eps=1e-4):
    left = df(lambda x: integral(f, 0, x), x)
    right = f(x)
    return abs(left - right) < eps
import math

def f(t):
    return t**2

print(theorem1(f, 2))   # 幾乎一定是 True
