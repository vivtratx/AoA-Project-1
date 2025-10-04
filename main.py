import csv
from collections import defaultdict

# ---------- geometry helpers ----------
def orient2D(a, b, c):
    # cross((b-a), (c-a)) ; >0: CCW, <0: CW
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

def signedArea(poly):
    s = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i+1) % len(poly)]
        s += x1*y2 - x2*y1
    return s

def ensureCCW(poly):
    if len(poly) >= 3 and signedArea(poly) < 0:
        poly.reverse()
    return poly

def nextIdx(i, n):  return (i + 1) % n
def prevIdx(i, n):  return (i - 1 + n) % n

def rightmostIndex(h):
    # max x, tie-break by smaller y
    return max(range(len(h)), key=lambda k: (h[k][0], -h[k][1]))

def leftmostIndex(h):
    # min x, tie-break by smaller y
    return min(range(len(h)), key=lambda k: (h[k][0],  h[k][1]))

# ---------- tiny base hull for n <= 3 (returns CCW) ----------
def smallHull(q):
    n = len(q)
    if n <= 2:
        return q[:]
    a, b, c = q
    s = orient2D(a, b, c)
    if s > 0:  return [a, b, c]      # CCW
    if s < 0:  return [a, c, b]      # flip to CCW
    # (no three collinear per spec, so this case won't occur)
    return [a, b, c]

# ---------- upper / lower tangents (strict tests; spec says no collinear) ----------
def upperTangent(L, R):
    """
    Return (i, j) for the upper tangent.
    Assumes L and R are CCW hulls, with L to the left of R.
    """
    nL, nR = len(L), len(R)
    i = rightmostIndex(L)
    j = leftmostIndex(R)
    while True:
        moved = False
        # move j forward on R while L[i] -> R[j] -> R[j+1] turns RIGHT (orient < 0)
        while orient2D(L[i], R[j], R[nextIdx(j, nR)]) < 0:
            j = nextIdx(j, nR); moved = True
        # move i backward on L while R[j] -> L[i] -> L[i-1] turns LEFT (orient > 0)
        while orient2D(R[j], L[i], L[prevIdx(i, nL)]) > 0:
            i = prevIdx(i, nL); moved = True
        if not moved:
            break
    return i, j

def lowerTangent(L, R):
    """
    Return (i, j) for the lower tangent.
    Assumes L and R are CCW hulls, with L to the left of R.
    """
    nL, nR = len(L), len(R)
    i = rightmostIndex(L)
    j = leftmostIndex(R)
    while True:
        moved = False
        # move j backward on R while L[i] -> R[j] -> R[j-1] turns LEFT (orient > 0)
        while orient2D(L[i], R[j], R[prevIdx(j, nR)]) > 0:
            j = prevIdx(j, nR); moved = True
        # move i forward on L while R[j] -> L[i] -> L[i+1] turns RIGHT (orient < 0)
        while orient2D(R[j], L[i], L[nextIdx(i, nL)]) < 0:
            i = nextIdx(i, nL); moved = True
        if not moved:
            break
    return i, j


# ---------- merge using tangents: walk *outer* arcs in CCW ----------
def mergeHullsUsingTangents(L, R, iU, jU, iL, jL):
    nL, nR = len(L), len(R)
    merged = []

    # L: CCW from lower -> upper (inclusive)
    i = iL
    merged.append(L[i])
    while i != iU:
        i = nextIdx(i, nL)
        merged.append(L[i])

    # R: CCW from upper -> lower (inclusive)
    j = jU
    merged.append(R[j])
    while j != jL:
        j = nextIdx(j, nR)
        merged.append(R[j])

    return merged  # CCW, no repeated start

# ---------- D&C hull on an x-sorted block ----------
def dcHull(sortedPts):
    n = len(sortedPts)
    if n <= 3:
        return ensureCCW(smallHull(sortedPts))
    mid = n // 2
    left  = ensureCCW(dcHull(sortedPts[:mid]))
    right = ensureCCW(dcHull(sortedPts[mid:]))
    iU, jU = upperTangent(left, right)
    iL, jL = lowerTangent(left, right)
    merged = mergeHullsUsingTangents(left, right, iU, jU, iL, jL)
    return ensureCCW(merged)

# ---------- public: hull indices with deterministic anchor ----------
def convexHullIndices(points):
    # map (x,y) -> original indices (handles duplicates, though spec says x are distinct)
    where = defaultdict(list)
    for idx, (x, y) in enumerate(points):
        where[(x, y)].append(idx)
    for k in where:
        where[k].sort(reverse=True)   # pop() is O(1)

    # Points are already x-sorted per spec; sorting again is harmless
    ptsSorted = sorted(points, key=lambda p: (p[0], p[1]))
    mid = len(ptsSorted) // 2

    # build top-level left/right hulls and tangents
    leftHull  = ensureCCW(dcHull(ptsSorted[:mid]))
    rightHull = ensureCCW(dcHull(ptsSorted[mid:]))
    iU, jU = upperTangent(leftHull, rightHull)
    iL, jL = lowerTangent(leftHull, rightHull)

    # deterministic anchor = left endpoint of top-level lower tangent
    anchorPoint = leftHull[iL]

    # final merged hull (top-level)
    hullPts = mergeHullsUsingTangents(leftHull, rightHull, iU, jU, iL, jL)

    # rotate so list starts at anchorPoint
    k = hullPts.index(anchorPoint)  # anchor must be on merged hull
    hullPts = hullPts[k:] + hullPts[:k]

    # map to original indices (0-based)
    hullIdx = [ where[(x, y)].pop() for (x, y) in hullPts ]
    return hullIdx

# ---------- I/O ----------
def readPointsCsv(path="input.csv"):
    pts = []
    with open(path, newline="") as f:
        rdr = csv.reader(f)
        for row in rdr:
            if not row: continue
            try:
                x = float(row[0].strip()); y = float(row[1].strip())
            except ValueError:
                continue   # skip header/bad rows
            pts.append([x, y])
    return pts

def writeOnePerLine(path, indices):
    with open(path, "w") as f:
        for i in indices:
            f.write(f"{i}\n")

# ---------- main ----------
if __name__ == "__main__":
    import sys
    inPath  = sys.argv[1] if len(sys.argv) > 1 else "input.csv"
    outPath = sys.argv[2] if len(sys.argv) > 2 else "output.txt"

    points = readPointsCsv(inPath)
    hullIdx = convexHullIndices(points)

    print("\n".join(str(i) for i in hullIdx))
    writeOnePerLine(outPath, hullIdx)
