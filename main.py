import csv
import sys
from typing import List, Tuple

Point = Tuple[float, float, int]                            # Represents points as (x, y, original_index)

def orient(a: Point, b: Point, c: Point) -> float:          # cross( b-a, c-a )
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

def is_ccw(a: Point, b: Point, c: Point) -> bool:           #checks if points are counter clockwise, true if a->b->c is a left turn
    return orient(a, b, c) > 0

def ccw_sort3(pts: List[Point]) -> List[Point]:             #sorts points in counter clockwise order, no collinear triples per spec
    a, b, c = pts
    if is_ccw(a, b, c): return [a, b, c]
    if is_ccw(a, c, b): return [a, c, b]
    if is_ccw(b, a, c): return [b, a, c]
    if is_ccw(b, c, a): return [b, c, a]
    if is_ccw(c, a, b): return [c, a, b]
    return [c, b, a]

def next_i(i: int, n: int) -> int: return (i + 1) % n        #move one step forward
def prev_i(i: int, n: int) -> int: return (i - 1) % n        #move one step backward

def rightmost_idx(poly: List[Point]) -> int:                 #checks for the rightmost index for d&c
    # max x, tie-breaker: higher y
    k = 0
    for i in range(1, len(poly)):
        if poly[i][0] > poly[k][0] or (poly[i][0] == poly[k][0] and poly[i][1] > poly[k][1]):
            k = i
    return k

def leftmost_idx(poly: List[Point]) -> int:                  #checks for the leftmost index for d&c
    # min x, tie-breaker: lower y
    k = 0
    for i in range(1, len(poly)):
        if poly[i][0] < poly[k][0] or (poly[i][0] == poly[k][0] and poly[i][1] < poly[k][1]):
            k = i
    return k

def getUpperTangent(L: List[Point], R: List[Point]) -> Tuple[int, int]:    #computes upper tangent between L and R
    i = rightmost_idx(L)
    j = leftmost_idx(R)
    nL, nR = len(L), len(R)
    while True:
        moved = False
        # advance i CCW while L[i+1] is above line (L[i], R[j])
        while orient(L[i], R[j], L[next_i(i, nL)]) > 0:
            i = next_i(i, nL); moved = True
        # advance j CW while R[j-1] is above line (L[i], R[j])
        while orient(R[j], L[i], R[prev_i(j, nR)]) < 0:
            j = prev_i(j, nR); moved = True
        if not moved:
            break
    return i, j

def getLowerTangent(L: List[Point], R: List[Point]) -> Tuple[int, int]:   #computes lower tanger between L and R
    i = rightmost_idx(L)
    j = leftmost_idx(R)
    nL, nR = len(L), len(R)
    while True:
        moved = False
        # advance i CW while L[i-1] is below line (L[i], R[j])
        while orient(L[i], R[j], L[prev_i(i, nL)]) < 0:
            i = prev_i(i, nL); moved = True
        # advance j CCW while R[j+1] is below line (L[i], R[j])
        while orient(R[j], L[i], R[next_i(j, nR)]) > 0:
            j = next_i(j, nR); moved = True
        if not moved:
            break
    return i, j

def merge_hulls(L: List[Point], R: List[Point]) -> List[Point]:          # Merges CCW hull L and CCW hull R into CCW hull
    if not L: return R[:]
    if not R: return L[:]

    iU, jU = getUpperTangent(L, R)
    iL, jL = getLowerTangent(L, R)

    merged = []

    # walk L from upper to lower CCW
    i = iU
    merged.append(L[i])
    while i != iL:
        i = next_i(i, len(L))
        merged.append(L[i])

    # walk R from lower to upper CCW
    j = jL
    merged.append(R[j])
    while j != jU:
        j = next_i(j, len(R))
        merged.append(R[j])

    return merged

def dac(points_sorted_by_x: List[Point]) -> List[Point]:  # Recursively compute the convex hull of points_sorted_by_x / doing basic divide and conqure
    n = len(points_sorted_by_x)
    if n <= 1: return points_sorted_by_x[:]
    if n == 2: return points_sorted_by_x[:]
    if n == 3: return ccw_sort3(points_sorted_by_x[:])
    mid = n // 2
    L = dac(points_sorted_by_x[:mid])
    R = dac(points_sorted_by_x[mid:])
    return merge_hulls(L, R)

def read_points_csv(path: str) -> List[Point]:           # Reads the points from the CSV file
    pts = []
    with open(path, "r") as f:
        for i, row in enumerate(csv.reader(f)):
            x, y = float(row[0]), float(row[1])
            pts.append((x, y, i))   # keep original index
    # input is supposed to be sorted by x; sort again as a guard
    pts.sort(key=lambda p: (p[0], p[1]))
    return pts

def write_indices(hull: List[Point], path: str) -> None:     # Writes the points to the output.txt file
    with open(path, "w") as f:
        for _, _, idx in hull:
            f.write(f"{idx}\n")

def main():                                                 #...It's main...idk what to tell you...
    in_path = sys.argv[1] if len(sys.argv) > 1 else "input.csv"   
    pts = read_points_csv(in_path)
    hull = dac(pts)
    write_indices(hull, "output.txt")

if __name__ == "__main__":
    main()
