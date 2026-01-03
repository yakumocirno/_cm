"""
幾何物件：點 Point、直線 Line、圓 Circle、三角形 Triangle
功能：
1) 兩直線交點
2) 兩圓交點
3) 直線與圓交點
4) 給定直線 L 與線外點 P，作垂線；求垂足 H
5) 用 (線上點 A、線外點 P、垂足 H) 形成的直角三角形驗證畢氏定理
6) 對上述物件做 平移/縮放/旋轉（縮放採「等比例縮放」，否則圓會變橢圓）

背後數學（都用註解寫在程式裡）：
- 向量、內積、投影：垂足與垂線
- 解析幾何：直線一般式 ax+by+c=0
- 聯立方程：直線-直線、直線-圓、圓-圓
- 旋轉矩陣、縮放、平移：仿射變換
"""

from __future__ import annotations
from dataclasses import dataclass
import math
from typing import List, Optional, Tuple

EPS = 1e-9

def is_close(a: float, b: float, eps: float = EPS) -> bool:
    return abs(a - b) <= eps * max(1.0, abs(a), abs(b))

def clamp0(x: float) -> float:
    # 讓 -1e-12 這種浮點誤差回到 0
    return 0.0 if abs(x) < 1e-12 else x

@dataclass(frozen=True)
class Point:
    x: float
    y: float

    # --- 向量運算 ---
    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, k: float) -> "Point":
        return Point(self.x * k, self.y * k)

    __rmul__ = __mul__

    def dot(self, other: "Point") -> float:
        # 內積：u·v = |u||v|cosθ，用來判斷垂直/投影
        return self.x * other.x + self.y * other.y

    def norm2(self) -> float:
        return self.x * self.x + self.y * self.y

    def dist2(self, other: "Point") -> float:
        return (self - other).norm2()

    def dist(self, other: "Point") -> float:
        return math.hypot(self.x - other.x, self.y - other.y)

    # --- 幾何變換（仿射） ---
    def translate(self, dx: float, dy: float) -> "Point":
        # 平移： (x, y) -> (x+dx, y+dy)
        return Point(self.x + dx, self.y + dy)

    def scale(self, s: float, about: "Point" = None) -> "Point":
        # 等比例縮放（以 about 為中心）：
        # p' = about + s*(p - about)
        if about is None:
            about = Point(0.0, 0.0)
        v = self - about
        return about + v * s

    def rotate(self, theta: float, about: "Point" = None) -> "Point":
        # 旋轉（以 about 為中心）：
        # [x']   [cos -sin][x]
        # [y'] = [sin  cos][y]
        if about is None:
            about = Point(0.0, 0.0)
        v = self - about
        c, s = math.cos(theta), math.sin(theta)
        return Point(about.x + c * v.x - s * v.y,
                     about.y + s * v.x + c * v.y)

@dataclass(frozen=True)
class Line:
    """
    直線一般式：a x + b y + c = 0
    法向量 n = (a, b)，方向向量 d = (b, -a)
    """
    a: float
    b: float
    c: float

    @staticmethod
    def from_points(p1: Point, p2: Point) -> "Line":
        # 兩點式轉一般式：
        # 通過 p1(x1,y1), p2(x2,y2)
        # 可取 a = y1 - y2, b = x2 - x1, c = x1*y2 - x2*y1
        a = p1.y - p2.y
        b = p2.x - p1.x
        c = p1.x * p2.y - p2.x * p1.y
        return Line(a, b, c).normalized()

    def normalized(self) -> "Line":
        # 讓 (a,b) 的長度為 1：利於距離/投影運算更穩定
        norm = math.hypot(self.a, self.b)
        if norm < EPS:
            raise ValueError("Invalid line: a and b are both ~0")
        return Line(self.a / norm, self.b / norm, self.c / norm)

    def direction(self) -> Point:
        # 方向向量 d = (b, -a)（與法向量垂直）
        return Point(self.b, -self.a)

    def eval(self, p: Point) -> float:
        return self.a * p.x + self.b * p.y + self.c

    def distance(self, p: Point) -> float:
        # 若已 normalized，距離 = |ax+by+c|
        # 一般情況距離 = |ax+by+c| / sqrt(a^2+b^2)
        denom = math.hypot(self.a, self.b)
        return abs(self.eval(p)) / denom

    def project_foot(self, p: Point) -> Point:
        """
        垂足 H：把點 p 投影到直線上
        直線 ax+by+c=0，法向量 n=(a,b)
        投影公式（解析幾何常用結論）：
            H = p - n * (ax_p + by_p + c) / (a^2 + b^2)
        """
        denom = self.a * self.a + self.b * self.b
        t = self.eval(p) / denom
        return Point(p.x - self.a * t, p.y - self.b * t)

    def perpendicular_through(self, p: Point) -> "Line":
        """
        經過點 p 且垂直於本線的直線：
        本線方向向量 d=(b,-a)
        垂直線的法向量可取 d（因為直線一般式的 (a,b) 是法向量）
        所以新線： (b) x + (-a) y + c2 = 0，代入 p 求 c2
        """
        A = self.b
        B = -self.a
        C = -(A * p.x + B * p.y)
        return Line(A, B, C).normalized()

    def two_points(self) -> Tuple[Point, Point]:
        """
        取直線上兩個不同點（方便做變換：把兩點變換後重建直線）。
        """
        if abs(self.b) > EPS:
            # 令 x=0, x=1 求 y
            y0 = -(self.c + self.a * 0.0) / self.b
            y1 = -(self.c + self.a * 1.0) / self.b
            return Point(0.0, y0), Point(1.0, y1)
        else:
            # b≈0 => ax + c = 0 => x = -c/a，取 y=0, y=1
            if abs(self.a) < EPS:
                raise ValueError("Invalid line")
            x = -self.c / self.a
            return Point(x, 0.0), Point(x, 1.0)

    # --- 交點 ---
    def intersect_line(self, other: "Line") -> Optional[Point]:
        """
        兩直線交點：解聯立
            a1 x + b1 y + c1 = 0
            a2 x + b2 y + c2 = 0
        用克拉瑪法則：
            D = a1*b2 - a2*b1
            x = (b1*c2 - b2*c1)/D
            y = (c1*a2 - c2*a1)/D
        平行或重合 => D≈0 => 無唯一交點
        """
        a1, b1, c1 = self.a, self.b, self.c
        a2, b2, c2 = other.a, other.b, other.c
        D = a1 * b2 - a2 * b1
        if abs(D) < EPS:
            return None
        x = (b1 * c2 - b2 * c1) / D
        y = (c1 * a2 - c2 * a1) / D
        return Point(x, y)

    # --- 幾何變換 ---
    def translate(self, dx: float, dy: float) -> "Line":
        # 平移：把線上兩點平移後再重建直線
        p1, p2 = self.two_points()
        return Line.from_points(p1.translate(dx, dy), p2.translate(dx, dy))

    def scale(self, s: float, about: Point = None) -> "Line":
        # 等比例縮放：同上（兩點法）
        p1, p2 = self.two_points()
        return Line.from_points(p1.scale(s, about), p2.scale(s, about))

    def rotate(self, theta: float, about: Point = None) -> "Line":
        p1, p2 = self.two_points()
        return Line.from_points(p1.rotate(theta, about), p2.rotate(theta, about))

@dataclass(frozen=True)
class Circle:
    center: Point
    r: float

    def __post_init__(self):
        if self.r < 0:
            raise ValueError("Radius must be non-negative")

    # --- 交點 ---
    def intersect_circle(self, other: "Circle") -> List[Point]:
        """
        兩圓交點（標準幾何推導）：
        圓心距 d = |C2-C1|
        若 d > r1+r2 或 d < |r1-r2| 或 (d≈0 且 r1≈r2) => 無(或無限多)交點
        否則在 C1->C2 方向上找到交弦中點，再用垂直方向偏移求兩交點。
        """
        c1, r1 = self.center, self.r
        c2, r2 = other.center, other.r
        dx, dy = c2.x - c1.x, c2.y - c1.y
        d = math.hypot(dx, dy)

        # 特殊情況
        if d < EPS and abs(r1 - r2) < EPS:
            return []  # 同圓：無限多交點，這裡回空表表示「不唯一」
        if d > r1 + r2 + EPS:
            return []
        if d < abs(r1 - r2) - EPS:
            return []
        if d < EPS:
            return []

        # a = 從 c1 沿著 c1->c2 到交弦中點的距離
        a = (r1*r1 - r2*r2 + d*d) / (2*d)
        # h = 從交弦中點到交點的垂直距離
        h2 = r1*r1 - a*a
        h2 = clamp0(h2)
        h = math.sqrt(h2)

        ux, uy = dx / d, dy / d  # 單位向量
        px = c1.x + a * ux
        py = c1.y + a * uy        # 交弦中點 P

        # 垂直方向 ( -uy, ux )
        rx = -uy * h
        ry =  ux * h

        pA = Point(px + rx, py + ry)
        pB = Point(px - rx, py - ry)

        if pA.dist2(pB) < 1e-16:
            return [pA]  # 相切
        return [pA, pB]

    def intersect_line(self, line: Line) -> List[Point]:
        """
        直線-圓交點：
        用「垂足+半弦長」最乾淨：
        - H = 圓心到直線的垂足（投影）
        - dist = 圓心到直線距離
        若 dist > r => 無交點
        若 dist = r => 相切，一點 H
        若 dist < r => 兩交點在直線方向上：H ± t * dir_unit
          t = sqrt(r^2 - dist^2)
        """
        L = line.normalized()
        O = self.center
        H = L.project_foot(O)
        dist = O.dist(H)

        if dist > self.r + EPS:
            return []
        if is_close(dist, self.r):
            return [H]

        t2 = self.r * self.r - dist * dist
        t2 = clamp0(t2)
        t = math.sqrt(t2)

        d = L.direction()
        dn = math.hypot(d.x, d.y)
        diru = Point(d.x / dn, d.y / dn)

        p1 = H + diru * t
        p2 = H - diru * t
        if p1.dist2(p2) < 1e-16:
            return [p1]
        return [p1, p2]

    # --- 幾何變換 ---
    def translate(self, dx: float, dy: float) -> "Circle":
        return Circle(self.center.translate(dx, dy), self.r)

    def scale(self, s: float, about: Point = None) -> "Circle":
        # 等比例縮放：半徑也乘 |s|
        return Circle(self.center.scale(s, about), abs(s) * self.r)

    def rotate(self, theta: float, about: Point = None) -> "Circle":
        # 旋轉不改變半徑
        return Circle(self.center.rotate(theta, about), self.r)

@dataclass(frozen=True)
class Triangle:
    A: Point
    B: Point
    C: Point

    def side_lengths(self) -> Tuple[float, float, float]:
        return (self.A.dist(self.B), self.B.dist(self.C), self.C.dist(self.A))

    def area_signed2(self) -> float:
        # 2倍有向面積：cross(B-A, C-A)
        return (self.B.x - self.A.x) * (self.C.y - self.A.y) - (self.B.y - self.A.y) * (self.C.x - self.A.x)

    def translate(self, dx: float, dy: float) -> "Triangle":
        return Triangle(self.A.translate(dx, dy), self.B.translate(dx, dy), self.C.translate(dx, dy))

    def scale(self, s: float, about: Point = None) -> "Triangle":
        return Triangle(self.A.scale(s, about), self.B.scale(s, about), self.C.scale(s, about))

    def rotate(self, theta: float, about: Point = None) -> "Triangle":
        return Triangle(self.A.rotate(theta, about), self.B.rotate(theta, about), self.C.rotate(theta, about))

# ---------------------------
# Demo：交點 + 垂線 + 畢氏定理驗證 + 變換示範
# ---------------------------
if __name__ == "__main__":
    print("=== 1) 兩直線交點 ===")
    L1 = Line.from_points(Point(0, 0), Point(2, 2))        # y=x
    L2 = Line.from_points(Point(0, 2), Point(2, 0))        # y=-x+2
    P = L1.intersect_line(L2)
    print("L1 ∩ L2 =", P)

    print("\n=== 2) 兩圓交點 ===")
    C1 = Circle(Point(0, 0), 2.0)
    C2 = Circle(Point(2, 0), 2.0)
    pts_cc = C1.intersect_circle(C2)
    print("C1 ∩ C2 =", pts_cc)

    print("\n=== 3) 直線與圓交點 ===")
    L3 = Line.from_points(Point(-3, 1), Point(3, 1))       # y=1
    pts_lc = C1.intersect_line(L3)
    print("C1 ∩ L3 =", pts_lc)

    print("\n=== 4) 給直線與線外點作垂線，求垂足 ===")
    L = Line.from_points(Point(-2, 0), Point(3, 0))        # x軸 y=0
    P0 = Point(1, 3)                                       # 線外點
    H = L.project_foot(P0)                                 # 垂足
    Lperp = L.perpendicular_through(P0)                    # 垂線
    print("Line L:", L)
    print("Point P0:", P0)
    print("Foot H:", H)
    print("Perpendicular line through P0:", Lperp)
    # 檢查：H 在 L 上、P0 在垂線上、且 (P0-H) ⟂ L 的方向
    print("Check H on L:", is_close(L.eval(H), 0.0))
    print("Check P0 on Lperp:", is_close(Lperp.eval(P0), 0.0))

    print("\n=== 5) 用直角三角形驗證畢氏定理 ===")
    # 取線上另一點 A（在 x軸上）
    A = Point(-4, 0)
    # 三角形頂點：P0, H, A，其中 ∠PHA 為直角（PH ⟂ HA）
    tri = Triangle(P0, H, A)
    PH2 = P0.dist2(H)
    HA2 = H.dist2(A)
    PA2 = P0.dist2(A)
    print("P0:", P0, "H:", H, "A:", A)
    print("PH^2 =", PH2)
    print("HA^2 =", HA2)
    print("PA^2 =", PA2)
    print("PH^2 + HA^2 =", PH2 + HA2)
    print("Pythagorean holds:", is_close(PH2 + HA2, PA2, eps=1e-7))

    print("\n=== 6) 幾何物件：平移/縮放/旋轉 ===")
    theta = math.radians(30)
    about = Point(0, 0)

    # Point
    p = Point(1, 2)
    print("Point p:", p)
    print("  translate:", p.translate(3, -1))
    print("  scale:", p.scale(2, about))
    print("  rotate:", p.rotate(theta, about))

    # Line
    print("Line L1:", L1)
    print("  translate:", L1.translate(1, 1))
    print("  scale:", L1.scale(1.5, about))
    print("  rotate:", L1.rotate(theta, about))

    # Circle
    print("Circle C1:", C1)
    print("  translate:", C1.translate(1, 1))
    print("  scale:", C1.scale(0.5, about))
    print("  rotate:", C1.rotate(theta, about))

    # Triangle
    print("Triangle tri:", tri)
    print("  translate:", tri.translate(1, 1))
    print("  scale:", tri.scale(2, about))
    print("  rotate:", tri.rotate(theta, about))
