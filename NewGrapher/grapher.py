import sys
import numpy as np 
import matplotlib.pyplot
import random

def readInPoints(file_name):
	'''
	Points are seperated by line, with x and y seperated by ','.
	'''
	f = open(file_name, "r")
	lines = f.readlines()
	points = np.zeros((3, len(lines)))
	i = 0
	for l in lines:
		l = l.split(",")
		points[:, i] = int(l[0].strip()), float(l[1].strip()), float(l[2].strip())
		i += 1
	return points

if __name__=="__main__":
	if (len(sys.argv) < 1):
		print("FORMAT: python", sys.argv[0], "<points file> <optional: output name>")
		print("For file formatting details, see source code for readIn functions.")
		exit(-1)
	fname = None if len(sys.argv) < 2 else sys.argv[2]
	# plot points
	points = readInPoints(sys.argv[1])

	matplotlib.pyplot.scatter(points[1], points[2], color='blue')

	ax = matplotlib.pyplot.gca()
	ax.set_aspect(0.4)
	
	for i in range(int(np.max(points[0]) + 1)):
		colors = ['red', 'blue', 'orange', 'purple']
		matplotlib.pyplot.scatter(points[1][points[0] == i], points[2][points[0] == i], color=colors[i])
		matplotlib.pyplot.plot(points[1][points[0] == i], points[2][points[0] == i], linestyle='dashed', color=colors[i])
	if fname is not None:
		matplotlib.pyplot.savefig(fname)
	matplotlib.pyplot.show()
