import sys
import numpy as np 
import matplotlib.pyplot


def readInPoints(file_name):
    '''
    Points are seperated by line, with x and y seperated by ','.
    '''
    f = open(file_name, "r")
    lines = f.readlines()
    points = np.zeros((2, len(lines)))
    i = 0
    for l in lines:
        l = l.split(",")
        points[:, i] = float(l[0].strip()), float(l[1].strip())
        i += 1
    return points

def readInHull(file_name):
    '''
    Assumes hull is constructed from lines built from sequential points
    in the file, then a final line wrapping around, connecting the last point to the first.
    '''
    f = open(file_name, "r")
    lines = f.readlines()
    points = np.zeros((2, len(lines) + 1))
    i = 0
    for l in lines:
        l = l.split(",")
        points[:, i] = float(l[0].strip()), float(l[1].strip())
        i += 1
    points[:, i] = points[:, 0]
    return points


def cross(a, b, c):
    return (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])


if __name__=="__main__":
    if (len(sys.argv) < 2):
        print("FORMAT: python", sys.argv[0], "<points file>")
        print("For file formatting details, see source code for readIn functions.")
        exit(-1)
    # plot points
    points = readInPoints(sys.argv[1])
    points = points[:,points[0,:].argsort()]
    numPoints = points.shape[1]
    bottomHull = []
    for p in range(numPoints):
        while len(bottomHull) >= 2 and cross(bottomHull[-2], bottomHull[-1], points[:, p]) <= 0:
            bottomHull.pop(-1)
        bottomHull.append(points[:, p])
    '''
    p = numPoints - 1
    topHull = []
    while p >= 0:
        while len(topHull) >= 2 and cross(topHull[-2], topHull[-1], points[:, p]) <= 0:
            topHull.pop(-1)
        topHull.append(points[:, p])
        p -= 1
    '''
    p = 0
    topHull = []
    while p < numPoints:
        while len(topHull) >= 2 and cross(topHull[-2], topHull[-1], points[:, p]) >= 0:
            topHull.pop(-1)
        topHull.append(points[:, p])
        p += 1
    # TODO : SUCESS, TESTING APPEARS TO SHOW SOLID RESULTS HERE, SO LONG AS Y's AREN"T A PROBLEM
    # SCRATCH OUT A THING FOR THIS ON PAPER!!!!!
    topHull.reverse()
    topHull.pop()
    bottomHull.pop()
    for h in (topHull + bottomHull):
        print(h.tolist())