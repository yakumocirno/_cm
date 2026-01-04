# -*- coding: utf-8 -*-
"""
Finite Field GF(p) implementation (p prime)
- Additive group: (F, +)
- Multiplicative group: (F\{0}, *)
- Distributivity: a*(b+c)=a*b+a*c

This file includes:
- FieldElement with operator overloading
- FiniteField
- AddGroup / MulGroup classes
- group axioms checker
- field distributivity checker
- demo usage

Copy & run directly.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, List, Iterable, Any, Optional


# =========================
# Basic utilities
# =========================

def is_prime(p: int) -> bool:
    if p <= 1:
        return False
    if p <= 3:
        return True
    if p % 2 == 0:
        return False
    i = 3
    while i * i <= p:
        if p % i == 0:
            return False
        i += 2
    return True


def egcd(a: int, b: int):
    """Extended GCD: returns (g, x, y) s.t. ax + by = g = gcd(a,b)"""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b != 0:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def modinv(a: int, p: int) -> int:
    """Inverse of a mod p (p prime, a != 0 mod p)."""
    a %= p
    if a == 0:
        raise ZeroDivisionError("0 has no multiplicative inverse in a field.")
    g, x, _ = egcd(a, p)
    if g != 1:
        # In GF(p) with prime p and a != 0, gcd must be 1
        raise ZeroDivisionError(f"{a} has no inverse mod {p} (gcd={g}).")
    return x % p


# =========================
# Field Element
# =========================

@dataclass(frozen=True)
class FieldElement:
    """An element of GF(p)"""
    value: int
    field: "FiniteField"

    def __post_init__(self):
        object.__setattr__(self, "value", self.value % self.field.p)

    def __int__(self) -> int:
        return self.value

    def __repr__(self) -> str:
        return f"{self.value} (mod {self.field.p})"

    def _coerce(self, other: Any) -> "FieldElement":
        if isinstance(other, FieldElement):
            if other.field != self.field:
                raise TypeError("Cannot operate on elements from different fields.")
            return other
        if isinstance(other, int):
            return self.field(other)
        raise TypeError(f"Unsupported operand type: {type(other)}")

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, FieldElement):
            return self.field == other.field and self.value == other.value
        if isinstance(other, int):
            return self.value == (other % self.field.p)
        return False

    # ---- addition / subtraction
    def __add__(self, other: Any) -> "FieldElement":
        o = self._coerce(other)
        return self.field(self.value + o.value)

    def __radd__(self, other: Any) -> "FieldElement":
        return self.__add__(other)

    def __neg__(self) -> "FieldElement":
        return self.field(-self.value)

    def __sub__(self, other: Any) -> "FieldElement":
        o = self._coerce(other)
        return self.field(self.value - o.value)

    def __rsub__(self, other: Any) -> "FieldElement":
        # other - self
        o = self._coerce(other)
        return self.field(o.value - self.value)

    # ---- multiplication / division
    def __mul__(self, other: Any) -> "FieldElement":
        o = self._coerce(other)
        return self.field(self.value * o.value)

    def __rmul__(self, other: Any) -> "FieldElement":
        return self.__mul__(other)

    def inv(self) -> "FieldElement":
        if self.value == 0:
            raise ZeroDivisionError("0 has no multiplicative inverse.")
        return self.field(modinv(self.value, self.field.p))

    def __truediv__(self, other: Any) -> "FieldElement":
        o = self._coerce(other)
        return self * o.inv()

    def __rtruediv__(self, other: Any) -> "FieldElement":
        # other / self
        o = self._coerce(other)
        return o * self.inv()


# =========================
# Finite Field GF(p)
# =========================

class FiniteField:
    def __init__(self, p: int):
        if not isinstance(p, int):
            raise TypeError("p must be an int.")
        if not is_prime(p):
            raise ValueError("GF(p) requires p to be prime in this implementation.")
        self.p = p

    def __call__(self, x: int) -> FieldElement:
        if not isinstance(x, int):
            raise TypeError("Only int can be coerced into GF(p) elements.")
        return FieldElement(x % self.p, self)

    def zero(self) -> FieldElement:
        return self(0)

    def one(self) -> FieldElement:
        return self(1)

    def elements(self) -> List[FieldElement]:
        return [self(i) for i in range(self.p)]

    def nonzero_elements(self) -> List[FieldElement]:
        return [self(i) for i in range(1, self.p)]

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, FiniteField) and self.p == other.p

    def __repr__(self) -> str:
        return f"GF({self.p})"


# =========================
# Group wrappers (like field_rational.py idea)
# =========================

class FiniteFieldAddGroup:
    """(F, +) over GF(p)"""
    def __init__(self, field: FiniteField):
        self.field = field

    def elements(self) -> List[FieldElement]:
        return self.field.elements()

    def op(self, a: FieldElement, b: FieldElement) -> FieldElement:
        return a + b

    def identity(self) -> FieldElement:
        return self.field.zero()

    def inv(self, a: FieldElement) -> FieldElement:
        return -a

    def eq(self, a: FieldElement, b: FieldElement) -> bool:
        return a == b


class FiniteFieldMulGroup:
    """(F\\{0}, *) over GF(p)"""
    def __init__(self, field: FiniteField):
        self.field = field

    def elements(self) -> List[FieldElement]:
        return self.field.nonzero_elements()

    def op(self, a: FieldElement, b: FieldElement) -> FieldElement:
        return a * b

    def identity(self) -> FieldElement:
        return self.field.one()

    def inv(self, a: FieldElement) -> FieldElement:
        return a.inv()

    def eq(self, a: FieldElement, b: FieldElement) -> bool:
        return a == b


# =========================
# group_axioms.py (checker)
# =========================

def check_group(
    G,
    *,
    check_commutative: bool = False,
    verbose: bool = True
) -> bool:
    """
    Check group axioms on a finite set:
    - closure
    - associativity
    - identity
    - inverse
    optionally commutative
    """
    elems = G.elements()
    e = G.identity()

    # closure
    for a in elems:
        for b in elems:
            c = G.op(a, b)
            if c not in elems:
                if verbose:
                    print("[FAIL] closure:", a, b, "->", c, "not in elements")
                return False

    # associativity
    for a in elems:
        for b in elems:
            for c in elems:
                left = G.op(G.op(a, b), c)
                right = G.op(a, G.op(b, c))
                if not G.eq(left, right):
                    if verbose:
                        print("[FAIL] associativity:", a, b, c, left, right)
                    return False

    # identity
    for a in elems:
        if not G.eq(G.op(a, e), a) or not G.eq(G.op(e, a), a):
            if verbose:
                print("[FAIL] identity:", a, "with e =", e)
            return False

    # inverse
    for a in elems:
        inva = G.inv(a)
        if not G.eq(G.op(a, inva), e) or not G.eq(G.op(inva, a), e):
            if verbose:
                print("[FAIL] inverse:", a, "inv =", inva, "e =", e)
            return False

    # commutative (abelian)
    if check_commutative:
        for a in elems:
            for b in elems:
                if not G.eq(G.op(a, b), G.op(b, a)):
                    if verbose:
                        print("[FAIL] commutative:", a, b)
                    return False

    if verbose:
        print("[OK] group axioms passed.")
    return True


# =========================
# field_axioms.py (checker)
# =========================

def check_distributivity(field: FiniteField, verbose: bool = True) -> bool:
    """
    Verify distributive law:
    a*(b+c) = a*b + a*c
    and (a+b)*c = a*c + b*c (both sides, though one implies other with commutativity of +)
    """
    elems = field.elements()

    for a in elems:
        for b in elems:
            for c in elems:
                left1 = a * (b + c)
                right1 = (a * b) + (a * c)
                if left1 != right1:
                    if verbose:
                        print("[FAIL] distributivity a*(b+c):", a, b, c, left1, right1)
                    return False

                left2 = (a + b) * c
                right2 = (a * c) + (b * c)
                if left2 != right2:
                    if verbose:
                        print("[FAIL] distributivity (a+b)*c:", a, b, c, left2, right2)
                    return False

    if verbose:
        print("[OK] distributivity passed.")
    return True


# =========================
# Demo
# =========================

def demo():
    F = FiniteField(7)          # GF(7)
    addG = FiniteFieldAddGroup(F)
    mulG = FiniteFieldMulGroup(F)

    print("Field:", F)
    print("Elements:", F.elements())
    print()

    # group checks
    print("Check additive group (should be abelian):")
    check_group(addG, check_commutative=True, verbose=True)
    print()

    print("Check multiplicative group (nonzero elements):")
    check_group(mulG, check_commutative=True, verbose=True)
    print()

    # field distributivity check
    print("Check distributivity:")
    check_distributivity(F, verbose=True)
    print()

    # operator overloading usage
    a = F(3)
    b = F(5)
    c = F(0)
    print("a =", a, "b =", b, "c =", c)
    print("a + b =", a + b)         # 3+5=8 mod7=1
    print("a - b =", a - b)         # 3-5=-2 mod7=5
    print("a * b =", a * b)         # 15 mod7=1
    print("b / a =", b / a)         # 5 * inv(3) mod7, inv(3)=5 since 3*5=15=1 mod7 -> 5*5=25=4
    print("b / a =", b / a)
    print("a == 3 ?", a == 3)
    print("a == 10 ?", a == 10)     # 10 mod7 = 3 -> True

    print("\nTry division by zero:")
    try:
        print(a / c)
    except ZeroDivisionError as e:
        print("Caught:", e)


if __name__ == "__main__":
    demo()
