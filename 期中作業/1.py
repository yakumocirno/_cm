def poly(c, z):
    s = 0
    for i, ci in enumerate(c):
        s += ci * z**i
    return s
