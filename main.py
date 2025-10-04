import csv
import matplotlib.pyplot as plt
import os

# START OF HELPERS
# Helper function to compute slope given 2 points in the form: p1 = [x, y] and p2 = [x, y]
def computeSlope(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    if x2-x1==0:
        raise ValueError("Undefined Slope")
    return (y2 - y1) / (x2 - x1)

def orient2D(a, b, c):
    # cross((b-a), (c-a)) ; >0: CCW, <0: CW, =0: collinear
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

def ensureCCW(poly):
    if len(poly) >= 3:
        area = 0.0
        for i in range(len(poly)):
            x1, y1 = poly[i]
            x2, y2 = poly[(i+1) % len(poly)]
            area += x1*y2 - x2*y1
        if area < 0:  # polygon is clockwise
            poly.reverse()
    return poly

def smallHull(q):
    n = len(q)
    if n <= 2:
        return q[:]                  # already sorted by x,y
    a, b, c = q
    o = orient2D(a, b, c)
    if o > 0:  return [a, b, c]      # CCW
    if o < 0:  return [a, c, b]      # make CCW
    return [a, c]                    # collinear -> keep endpoints only

def upperTangent(L, R):
    """
    Return (i, j): indices on L and R where the upper tangent touches.
    L and R are CCW convex hulls, points as [x, y].
    """
    nL, nR = len(L), len(R)

    # rightmost on L (max x, tie-break: smaller y)
    i = 0
    for k in range(1, nL):
        if (L[k][0] > L[i][0]) or (L[k][0] == L[i][0] and L[k][1] < L[i][1]):
            i = k

    # leftmost on R (min x, tie-break: smaller y)
    j = 0
    for k in range(1, nR):
        if (R[k][0] < R[j][0]) or (R[k][0] == R[j][0] and R[k][1] < R[j][1]):
            j = k

    # hill-climb to the upper tangent
    changed = True
    while changed:
        changed = False
        # move i CCW on L while turn (R[j] -> L[i] -> L[i+1]) is left
        while True:
            ip1 = (i + 1) % nL
            ax, ay = R[j]; bx, by = L[i]; cx, cy = L[ip1]
            cross = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
            if cross >= 0:
                i = ip1; changed = True
            else:
                break
        # move j CW on R while turn (L[i] -> R[j] -> R[j-1]) is right
        while True:
            jm1 = (j - 1) % nR
            ax, ay = L[i]; bx, by = R[j]; cx, cy = R[jm1]
            cross = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
            if cross <= 0:
                j = jm1; changed = True
            else:
                break
    return i, j

def lowerTangent(L, R):
    """
    Return (i, j): indices on L and R where the lower tangent touches.
    L and R are CCW convex hulls, points as [x, y].
    """
    nL, nR = len(L), len(R)

    # rightmost on L (max x, tie-break: smaller y)
    i = 0
    for k in range(1, nL):
        if (L[k][0] > L[i][0]) or (L[k][0] == L[i][0] and L[k][1] < L[i][1]):
            i = k

    # leftmost on R (min x, tie-break: smaller y)
    j = 0
    for k in range(1, nR):
        if (R[k][0] < R[j][0]) or (R[k][0] == R[j][0] and R[k][1] < R[j][1]):
            j = k

    # hill-climb to the lower tangent
    changed = True
    while changed:
        changed = False
        # move i CW on L while turn (R[j] -> L[i] -> L[i-1]) is right
        while True:
            im1 = (i - 1) % nL
            ax, ay = R[j]; bx, by = L[i]; cx, cy = L[im1]
            cross = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
            if cross <= 0:
                i = im1; changed = True
            else:
                break
        # move j CCW on R while turn (L[i] -> R[j] -> R[j+1]) is left
        while True:
            jp1 = (j + 1) % nR
            ax, ay = L[i]; bx, by = R[j]; cx, cy = R[jp1]
            cross = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
            if cross >= 0:
                j = jp1; changed = True
            else:
                break
    return i, j

# END OF HELPERS

# Functions that compute what we need
### Divide the set of points into two “halves” A and B according to x-coordinate.
def divide(points):
    mid = len(points)//2
    a = points[:mid]
    b = points[mid:]
    return a, b, points

### Recursively compute the convex hull of A, and re-cursively compute the convex hull of B.

# Let T be the line that connects A<->B where A is the rightmost point in CH(A) and B is the leftmost point in CH(B)
def dcHull(pts):
    n = len(pts)
    if n <= 3:
        return ensureCCW(smallHull(pts))
    mid = n // 2
    left  = ensureCCW(dcHull(pts[:mid]))
    right = ensureCCW(dcHull(pts[mid:]))

    merged = mergeHulls(left, right)   # expects CCW inputs
    return ensureCCW(merged)


### Merge the two convex hulls to get the convex hull of A ∪ B.
def mergeHulls(L, R):
    # Robust merge: take the union of hull vertices and run monotone chain
    pts = L + R

    # dedupe coordinates (keep one of each (x,y))
    seen = set()
    uniq = []
    for x, y in pts:
        if (x, y) in seen: 
            continue
        seen.add((x, y))
        uniq.append([x, y])

    # monotone chain on the small set 'uniq'
    uniq.sort(key=lambda p: (p[0], p[1]))

    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    lower = []
    for p in uniq:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(uniq):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    merged = lower[:-1] + upper[:-1]   # CCW hull
    return merged


### Merging: compute upper and lower tangent to the con-vex hulls in O(n) time
def convexHullIndices(points):
    # build a multi-map from (x,y) -> list of original indices (handles duplicates)
    from collections import defaultdict
    where = defaultdict(list)
    for idx, (x, y) in enumerate(points):
        where[(x, y)].append(idx)
    for key in where:
        where[key].sort(reverse=True)  # pop from end (O(1))

    # sort by x,y and run D&C
    ptsSorted = sorted(points, key=lambda p: (p[0], p[1]))
    hullPts = dcHull(ptsSorted)   # hull as [[x,y], ...] in CCW order

    # map hull points back to original indices
    result = []
    for x, y in hullPts:
        result.append(where[(x, y)].pop())
    return result

def rotateToStart(indices, startAt):
    if not indices: 
        return indices
    if startAt not in indices:
        return indices
    k = indices.index(startAt)
    return indices[k:] + indices[:k]

# Main Logic goes here:
points = []
x = []
y = []

with open('input.csv', newline='') as file:
    rdr = csv.reader(file)
    for row in rdr:
        if not row:
            continue
        try:
            x = float(row[0].strip())
            y = float(row[1].strip())
        except ValueError:
            continue
        points.append([x, y])

hullIdx = convexHullIndices(points)
hullIdx = rotateToStart(hullIdx, 22)
print("\n".join(str(i) for i in hullIdx))
