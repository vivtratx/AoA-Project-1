import csv
import matplotlib.pyplot as plt
import os


# Helper function to compute slope given 2 points in the form: p1 = [x, y] and p2 = [x, y]
def computeSlope(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    if x2-x1==0:
        raise ValueError("Undefined Slope")
    return (y2 - y1) / (x2 - x1)

def getUpperTangent(L, R):
    """
    Compute the upper tangent between two CCW convex hulls L and R.
    Returns (i, j): indices into L and R where the upper tangent touches.
    No inner helpers; everything is written inline for readability.
    """
    nL, nR = len(L), len(R)

    # Find rightmost point of L (max x, tie-breaker: smaller y)
    i = 0
    for k in range(1, nL):
        if (L[k][0] > L[i][0]) or (L[k][0] == L[i][0] and L[k][1] < L[i][1]):
            i = k

    # Find leftmost point of R (min x, tie-breaker: smaller y)
    j = 0
    for k in range(1, nR):
        if (R[k][0] < R[j][0]) or (R[k][0] == R[j][0] and R[k][1] < R[j][1]):
            j = k

    # Iterate until no endpoint moves
    changed = True
    while changed:
        changed = False

        # Move i CCW on L while turn (R[j] -> L[i] -> L[i+1]) is left (cross > 0)
        while True:
            ip1 = (i + 1) % nL
            ax, ay = R[j]
            bx, by = L[i]
            cx, cy = L[ip1]
            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross > 0:
                i = ip1
                changed = True
            else:
                break

        # Move j CW on R while turn (L[i] -> R[j] -> R[j-1]) is right (cross < 0)
        while True:
            jm1 = (j - 1) % nR
            ax, ay = L[i]
            bx, by = R[j]
            cx, cy = R[jm1]
            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross < 0:
                j = jm1
                changed = True
            else:
                break

    return i, j

def getLowerTangent(L, R):
    """
    Compute the lower tangent between two CCW convex hulls L and R.
    Returns (i, j): indices into L and R where the lower tangent touches.
    No inner helpers; everything is written inline for readability.
    """
    nL, nR = len(L), len(R)

    # Find rightmost point of L (max x, tie-breaker: smaller y)
    i = 0
    for k in range(1, nL):
        if (L[k][0] > L[i][0]) or (L[k][0] == L[i][0] and L[k][1] < L[i][1]):
            i = k

    # Find leftmost point of R (min x, tie-breaker: smaller y)
    j = 0
    for k in range(1, nR):
        if (R[k][0] < R[j][0]) or (R[k][0] == R[j][0] and R[k][1] < R[j][1]):
            j = k

    # Iterate until no endpoint moves
    changed = True
    while changed:
        changed = False

        # Move i CW on L while turn (R[j] -> L[i] -> L[i-1]) is right (cross < 0)
        while True:
            im1 = (i - 1) % nL
            ax, ay = R[j]
            bx, by = L[i]
            cx, cy = L[im1]
            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross < 0:
                i = im1
                changed = True
            else:
                break

        # Move j CCW on R while turn (L[i] -> R[j] -> R[j+1]) is left (cross > 0)
        while True:
            jp1 = (j + 1) % nR
            ax, ay = L[i]
            bx, by = R[j]
            cx, cy = R[jp1]
            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross > 0:
                j = jp1
                changed = True
            else:
                break

    return i, j

# Load all the points into an array of points
points = []
x = []
y = []

with open('input.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for line in csvFile:
        points.append(line)
        x.append(line[0])
        y.append(line[1])


# Main Logic goes down here:

### Divide the set of points into two “halves” A and B according to x-coordinate.

mid = len(points)//2
a = points[:mid]
b = points[mid:]

### Recursively compute the convex hull of A, and re-cursively compute the convex hull of B.

# Let T be the line that connects A<->B where A is the rightmost point in CH(A) and B is the leftmost point in CH(B)
while (T is not lower tangent to all):
    while(T is not lower tangent to A):
        Move a clockwise one position
        a = (a-1) % size(CHA(A))
        a--
    
    while(T is not lower tangent to B):
        Move b counterclockwise one position
        b = (b-1) % size(CHA(B))
        b++

### Merge the two convex hulls to get the convex hull of A ∪ B.
### Merging: compute upper and lower tangent to the con-vex hulls in O(n) time
