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

if __name__=="__main__":
	if (len(sys.argv) < 3):
		print("FORMAT: python", sys.argv[0], "<points file> <hull file>")
		print("For file formatting details, see source code for readIn functions.")
		exit(-1)
	# plot points
	points = readInPoints(sys.argv[1])
	matplotlib.pyplot.scatter(points[0], points[1]) 
	# plot hull
	hull = readInHull(sys.argv[2])
	matplotlib.pyplot.scatter(hull[0], hull[1], color='red')
	matplotlib.pyplot.show()