# -*- coding: utf-8 -*-
import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

EPS = 1e-9

def is_close(a: float, b: float, eps: float = EPS) -> bool:
    return abs(a - b) <= eps

# =========================
# Point
# =========================
@dataclass(frozen=True)
class Point:
    x: float
    y: float

    # --- vector ops ---
    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, k: float) -> "Point":
        return Point(self.x * k, self.y * k)

    def __rmul__(self, k: float) -> "Point":
        return self.__mul__(k)

    def dot(self, other: "Point") -> float:
        return self.x * other.x + self.y * other.y

    def cross(self, other: "Point") -> float:
        # 2D cross product magnitude (as scalar)
        return self.x * other.y - self.y * other.x

    def norm2(self) -> float:
        return self.x * self.x + self.y * self.y

    def norm(self) -> float:
        return math.sqrt(self.norm2())

    def dist(self, other: "Point") -> float:
        return (self - other).norm()

    # --- transforms ---
    def translate(self, dx: float, dy: float) -> "Point":
        return Point(self.x + dx, self.y + dy)

    def scale(self, s: float, center: "Point" = None) -> "Point":
        # scaling about center:  c + s*(p-c)
        if center is None:
            center = Point(0.0, 0.0)
        v = self - center
        return center + v * s

    def rotate(self, theta_rad: float, center: "Point" = None) -> "Point":
        # rotation about center:
        # p' = c + R(theta)*(p-c)
        if center is None:
            center = Point(0.0, 0.0)
        v = self - center
        c = math.cos(theta_rad)
        s = math.sin(theta_rad)
        v2 = Point(c * v.x - s * v.y, s * v.x + c * v.y)
        return center + v2

    def __repr__(self) -> str:
        return f"Point({self.x:.6f}, {self.y:.6f})"


# =========================
# Line: represented as ax + by + c = 0
# =========================
@dataclass(frozen=True)
class Line:
    a: float
    b: float
    c: float

    @staticmethod
    def from_points(p1: Point, p2: Point) -> "Line":
        # Through p1(x1,y1), p2(x2,y2)
        # a = y1 - y2, b = x2 - x1, c = x1*y2 - x2*y1
        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y
        a = y1 - y2
        b = x2 - x1
        c = x1 * y2 - x2 * y1
        return Line(a, b, c)

    def direction(self) -> Point:
        # direction vector along the line is (b, -a)
        return Point(self.b, -self.a)

    def normal(self) -> Point:
        return Point(self.a, self.b)

    def eval(self, p: Point) -> float:
        return self.a * p.x + self.b * p.y + self.c

    def is_parallel(self, other: "Line") -> bool:
        # a1*b2 - a2*b1 == 0
        return abs(self.a * other.b - other.a * self.b) <= EPS

    def intersection(self, other: "Line") -> Optional[Point]:
        # Solve:
        # a1 x + b1 y + c1 = 0
        # a2 x + b2 y + c2 = 0
        det = self.a * other.b - other.a * self.b
        if abs(det) <= EPS:
            return None  # parallel or coincident
        x = (self.b * other.c - other.b * self.c) / det
        y = (other.a * self.c - self.a * other.c) / det
        return Point(x, y)

    def projection_of(self, p: Point) -> Point:
        # Foot of perpendicular from p to this line:
        # line: ax+by+c=0
        # proj = ( x - a*(ax+by+c)/(a^2+b^2), y - b*(ax+by+c)/(a^2+b^2) )
        denom = self.a * self.a + self.b * self.b
        if denom <= EPS:
            raise ValueError("Invalid line: a=b=0")
        t = (self.a * p.x + self.b * p.y + self.c) / denom
        return Point(p.x - self.a * t, p.y - self.b * t)

    def perpendicular_through(self, p: Point) -> "Line":
        # A perpendicular line has normal vector parallel to this line's direction (b,-a)
        # Using normal n' = (b, -a):  b x + (-a) y + c' = 0, passing through p:
        a2 = self.b
        b2 = -self.a
        c2 = -(a2 * p.x + b2 * p.y)
        return Line(a2, b2, c2)

    # --- transforms ---
    def translate(self, dx: float, dy: float) -> "Line":
        # Substitute x = x' - dx, y = y' - dy:
        # a(x'-dx)+b(y'-dy)+c = 0 => a x' + b y' + (c - a dx - b dy)=0
        return Line(self.a, self.b, self.c - self.a * dx - self.b * dy)

    def scale(self, s: float, center: Point = None) -> "Line":
        # Transform two points on line then rebuild
        if center is None:
            center = Point(0.0, 0.0)
        p1, p2 = self.sample_points()
        p1s = p1.scale(s, center)
        p2s = p2.scale(s, center)
        return Line.from_points(p1s, p2s)

    def rotate(self, theta_rad: float, center: Point = None) -> "Line":
        if center is None:
            center = Point(0.0, 0.0)
        p1, p2 = self.sample_points()
        p1r = p1.rotate(theta_rad, center)
        p2r = p2.rotate(theta_rad, center)
        return Line.from_points(p1r, p2r)

    def sample_points(self) -> Tuple[Point, Point]:
        # Get two distinct points on the line for transform purposes
        # If b != 0, use x=0 and x=1
        if abs(self.b) > EPS:
            y0 = -(self.c) / self.b
            y1 = -(self.a * 1.0 + self.c) / self.b
            return Point(0.0, y0), Point(1.0, y1)
        # else b == 0 => a x + c=0 => x = -c/a
        if abs(self.a) > EPS:
            x0 = -(self.c) / self.a
            return Point(x0, 0.0), Point(x0, 1.0)
        raise ValueError("Invalid line: a=b=0")

    def __repr__(self) -> str:
        return f"Line({self.a:.6f}x + {self.b:.6f}y + {self.c:.6f} = 0)"


# =========================
# Circle
# =========================
@dataclass(frozen=True)
class Circle:
    center: Point
    r: float

    def contains(self, p: Point) -> bool:
        return abs(self.center.dist(p) - self.r) <= 1e-6

    # --- transforms ---
    def translate(self, dx: float, dy: float) -> "Circle":
        return Circle(self.center.translate(dx, dy), self.r)

    def scale(self, s: float, center: Point = None) -> "Circle":
        if center is None:
            center = Point(0.0, 0.0)
        return Circle(self.center.scale(s, center), abs(s) * self.r)

    def rotate(self, theta_rad: float, center: Point = None) -> "Circle":
        if center is None:
            center = Point(0.0, 0.0)
        return Circle(self.center.rotate(theta_rad, center), self.r)

    def __repr__(self) -> str:
        return f"Circle(center={self.center}, r={self.r:.6f})"


# =========================
# Triangle
# =========================
@dataclass(frozen=True)
class Triangle:
    a: Point
    b: Point
    c: Point

    def side_lengths(self) -> Tuple[float, float, float]:
        return (self.a.dist(self.b), self.b.dist(self.c), self.c.dist(self.a))

    def area(self) -> float:
        # 0.5 * |(b-a) x (c-a)|
        return 0.5 * abs((self.b - self.a).cross(self.c - self.a))

    # --- transforms ---
    def translate(self, dx: float, dy: float) -> "Triangle":
        return Triangle(self.a.translate(dx, dy),
                        self.b.translate(dx, dy),
                        self.c.translate(dx, dy))

    def scale(self, s: float, center: Point = None) -> "Triangle":
        if center is None:
            center = Point(0.0, 0.0)
        return Triangle(self.a.scale(s, center),
                        self.b.scale(s, center),
                        self.c.scale(s, center))

    def rotate(self, theta_rad: float, center: Point = None) -> "Triangle":
        if center is None:
            center = Point(0.0, 0.0)
        return Triangle(self.a.rotate(theta_rad, center),
                        self.b.rotate(theta_rad, center),
                        self.c.rotate(theta_rad, center))

    def __repr__(self) -> str:
        return f"Triangle({self.a}, {self.b}, {self.c})"


# =========================
# Intersections
# =========================
def intersect_line_circle(line: Line, circle: Circle) -> List[Point]:
    # Solve line + circle:
    # line: ax+by+c=0
    # circle: (x-h)^2 + (y-k)^2 = r^2
    # Use projection + distance method:
    h, k = circle.center.x, circle.center.y
    r = circle.r

    # distance from center to line:
    denom = math.sqrt(line.a * line.a + line.b * line.b)
    if denom <= EPS:
        raise ValueError("Invalid line for intersection.")
    d = abs(line.a * h + line.b * k + line.c) / denom

    if d > r + 1e-9:
        return []  # no intersection
    foot = line.projection_of(circle.center)

    if is_close(d, r, 1e-7):
        return [foot]  # tangent

    # Two points: move along line direction unit vector by t = sqrt(r^2 - d^2)
    dirv = line.direction()
    L = dirv.norm()
    if L <= EPS:
        raise ValueError("Invalid line direction.")
    u = Point(dirv.x / L, dirv.y / L)
    t = math.sqrt(max(0.0, r * r - d * d))
    p1 = foot + u * t
    p2 = foot - u * t
    return [p1, p2]


def intersect_circle_circle(c1: Circle, c2: Circle) -> List[Point]:
    # Standard circle-circle intersection geometry
    p0, r0 = c1.center, c1.r
    p1, r1 = c2.center, c2.r
    d = p0.dist(p1)

    # No solutions: separate, one inside other, or coincident
    if d > r0 + r1 + 1e-9:
        return []
    if d < abs(r0 - r1) - 1e-9:
        return []
    if d <= EPS and abs(r0 - r1) <= 1e-9:
        # Infinite intersections (same circle) -> return empty, ambiguous
        return []

    # a = distance from p0 to line of intersection points along p0->p1
    a = (r0 * r0 - r1 * r1 + d * d) / (2 * d)
    h2 = r0 * r0 - a * a
    if h2 < 0 and h2 > -1e-8:
        h2 = 0.0
    if h2 < 0:
        return []

    # point p2 is base point
    v = (p1 - p0) * (1.0 / d)  # unit vector from p0 to p1
    p2 = p0 + v * a

    if is_close(h2, 0.0, 1e-7):
        return [p2]  # tangent

    h = math.sqrt(h2)
    # perpendicular unit vector
    perp = Point(-v.y, v.x)
    i1 = p2 + perp * h
    i2 = p2 - perp * h
    return [i1, i2]


# =========================
# Pythagorean verification
# =========================
def verify_pythagoras(line: Line, point_on_line: Point, point_off_line: Point) -> Tuple[Triangle, float, float, float, bool]:
    # foot of perpendicular from off point to line
    foot = line.projection_of(point_off_line)
    tri = Triangle(point_on_line, foot, point_off_line)

    # lengths: let A=point_on_line, B=foot, C=off
    A, B, C = point_on_line, foot, point_off_line
    AB = A.dist(B)
    BC = B.dist(C)
    AC = A.dist(C)

    ok = is_close(AB * AB + BC * BC, AC * AC, 1e-6)
    return tri, AB, BC, AC, ok


# =========================
# Demo / Usage
# =========================
def main():
    print("=== 1) 定義點、線、圓 ===")
    P = Point(1, 2)
    Q = Point(5, 4)
    line1 = Line.from_points(Point(0, 0), Point(4, 4))   # y=x
    line2 = Line.from_points(Point(0, 4), Point(4, 0))   # y=-x+4
    circle1 = Circle(center=Point(0, 0), r=5)
    circle2 = Circle(center=Point(6, 0), r=5)

    print("P =", P)
    print("Q =", Q)
    print("line1 =", line1)
    print("line2 =", line2)
    print("circle1 =", circle1)
    print("circle2 =", circle2)

    print("\n=== 2) 兩直線交點 ===")
    I = line1.intersection(line2)
    print("intersection(line1, line2) =", I)  # should be (2,2)

    print("\n=== 3) 直線與圓交點 ===")
    pts_lc = intersect_line_circle(line1, circle1)  # y=x with radius 5
    print("intersect(line1, circle1) =", pts_lc)

    print("\n=== 4) 兩圓交點 ===")
    pts_cc = intersect_circle_circle(circle1, circle2)
    print("intersect(circle1, circle2) =", pts_cc)

    print("\n=== 5) 給定直線與線外一點，作垂線並求垂足 ===")
    off = Point(3, 0)              # not on y=x
    foot = line1.projection_of(off)
    perp = line1.perpendicular_through(off)
    print("off point =", off)
    print("foot on line1 =", foot)
    print("perpendicular line through off =", perp)
    # check perpendicular by dot(direction vectors)=0
    d1 = line1.direction()
    d2 = perp.direction()
    print("direction dot =", d1.dot(d2), "(~0 means perpendicular)")

    print("\n=== 6) 用 (線上點、垂足、線外點) 驗證畢氏定理 ===")
    on = Point(4, 4)               # on y=x
    tri, AB, BC, AC, ok = verify_pythagoras(line1, on, off)
    print("triangle =", tri)
    print("AB(on-foot) =", AB)
    print("BC(foot-off) =", BC)
    print("AC(on-off) =", AC)
    print("Check: AB^2 + BC^2 == AC^2 ? ->", ok)
    print("AB^2 + BC^2 =", AB*AB + BC*BC, " AC^2 =", AC*AC)

    print("\n=== 7) 定義三角形物件：邊長與面積 ===")
    t = Triangle(Point(0,0), Point(3,0), Point(0,4))
    print("t =", t)
    print("side lengths =", t.side_lengths())
    print("area =", t.area(), "(should be 6)")

    print("\n=== 8) 平移 / 縮放 / 旋轉 示範 ===")
    dx, dy = 2, -1
    theta = math.pi / 6  # 30 degrees
    center = Point(0, 0)

    print("P translate:", P.translate(dx, dy))
    print("P scale about origin (s=2):", P.scale(2, center))
    print("P rotate about origin (30deg):", P.rotate(theta, center))

    print("line1 translate:", line1.translate(dx, dy))
    print("line1 scale about origin:", line1.scale(2, center))
    print("line1 rotate about origin:", line1.rotate(theta, center))

    print("circle1 translate:", circle1.translate(dx, dy))
    print("circle1 scale about origin:", circle1.scale(2, center))
    print("circle1 rotate about origin:", circle1.rotate(theta, center))

    print("triangle t translate:", t.translate(dx, dy))
    print("triangle t scale:", t.scale(2, center))
    print("triangle t rotate:", t.rotate(theta, center))

if __name__ == "__main__":
    main()
