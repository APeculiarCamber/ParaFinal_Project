We have 3 algorithms implemented in parallel for computing the convex hull of a set of 2D points:
Chan's Algorithm in the /Chan_Algorithm directory
Quickhull in the /Quickhull directory
Monotone Chain in the /Monotone_Chain directory


POINTS.BIN:
All programs assume that there is a ./pts directory with a "points.bin" file, which is a binary file containing packed 2D points.
See Point Generation Use.


Each C program file for the algorithms is accompanied by 4 shell scripts to run the program:
"run_*.sh" compiles and runs the program to output the time taken to compute the convex hull, no output of the convex hull is given except its size.
"run_*_debug.sh" compiles and runs the program and outputs the time taken AND the convex hull.
"run_*_bulk.sh" runs the program assuming that it has already been compiled.
"run-all-*-bench.sh" runs the benchmarking tests 

"run-all-*-bench.sh" takes no command line arguments. 
On a per algorithm basis, the other 3 scripts take the same command line arguments, which are as follows for each algorithm:
-- Chan's Algorithm:
	run_Chan.sh <number of nodes> <number of ranks PER NODE> <TOTAL number of points>
-- Quickhull:
	run-quick.sh <number of nodes> <number of ranks PER NODE> <TOTAL number of points>
-- Monotone Chain:
	run_mono.sh <number of nodes> <number of ranks PER NODE> <TOTAL number of points> <oversampling>

-- OVERSAMPLING NOTES:
	-- oversampling is a special argument for monotone chain which specifies the amount of oversampling to perform for samplesort.
	-- oversampling should be a positive number which is less than (TOTAL number of points / TOTAL number of ranks).
	-- For our benchmarking case we used an oversampling of 64 for all cases. For our correctness testing on a data set of 2048, we used 8 or 4.


POINT GENERATION USE:
Point Generation files are contained within the /Points_Creator directory. The program can uniformly sample points from
an area bounded by a rectangle or by a circle.
The shell script "run_points.sh" compiles and runs the points generation program for generating the specified number of points 
and packs them into a binary file called "points.bin" within the working directory. 
The shell script "run_points_debug.sh" compiles and runs the points generation program for similarly generating points but also 
outputs the the points is a plain-text format to the slurm output of the job.
For the algorithm programs to run, it is required that there exists a pts/points.bin file so it is recommended that Points_Creator files
are copied to "ALGORITHM_HOME/pts" and run there, where ALGORITHM_HOME is the location you intend to run our algorithms from.
The command line arguments for both scripts is as follows:
	For circle bounding boxes:
	./run_points.sh <number of nodes> <number of ranks per NODE> <number of points PER RANK> <cuda thread count PER RANK> <radius>
	For rectangle bounding boxes:
	./run_points.sh <number of nodes> <number of ranks per NODE> <number of points PER RANK> <cuda thread count PER RANK> <left> <lower> <right> <upper>

COMMAND LINE NOTES:
	-- left < right, lower < upper. (left, lower) is the bottom left corner and (right, upper) the upper right corner.
	-- The threads per node shouldn't exceeed 256 per GPU (the argument is threads per rank). We allocate 4 GPUs.
	-- The number of ranks shouldn't exceed 32.
	-- VERY IMPORTANT: The number of points is PER RANK, so you will create (number of ranks) * (number of points) TOTAL POINTS.
		-- This is in contrast to our convex hull programs which expect the TOTAL number of points.
	-- For benchmarking, we used the command ./run_points.sh 2 32 16777216 8 0 0 1 1 for rectangle bounding, this creates 2^30 points total. Programs are able to read in fewer points than are available; we don't have to continuously re-generate points if we do this.



POINTS GRAPHER:
"grapher.py" is a python program which is used to plot a set of points and its convex hull and is found in /Grapher.
It generates a matplot with all points marked in blue and the convex hull marked in red. 
It requires matplot and numpy to function properly.
IT SHOULD ONLY BE RUN ON YOUR LOCAL MACHINE.
It's command line arguments are of the following format:
	-- python grapher.py <all points file> <hull points file> <optional: save plot to this file name>
Example run:
	-- python grapher.py example_pts.txt example_hull.txt
Notes:
	-- The input files are plain-text files with row seperated points, with their x and y components seperated by a comma. See the example text files.


VERIFYING CORRECTNESS RESULTS:
-- In the pts directory, run ./run_points_debug.sh for a managable number of points; we used 2048 for our correctness testing.
-- Open the slurm output for the job and copy the printed points to a text file on your local machine, remove any text that is not a points x or y value or the seperating comma. The grapher.py program only expects these 3 elements to exist on a line of an input file.
-- For any algorithm, run ./run_*_debug.sh, open the slurm output and copy the convex hull points printed to a text file on your local machine, similarly remove all elements from lines except the x component, the y component, and the separating comma.
-- run "python grapher.py <your points file> <your hull file>" to generate a plot of your points and the computed convex hull.


VERIFYING BENCHMARKING RESULTS:
-- Run the associated ./run-all-*-bench.sh script and wait for the jobs to finish.
-- grep the slurm outputs for the times taken using 'grep "Time " *'
-- In the majority of cases, grep will provide the times in the order the jobs were requested, see the ./run-all-*-bench.sh script.



If it's not overstepping, feel free to email idemaj@rpi.edu with questions about running programs. 
Use may be tricky but we are confident our algorithms work.