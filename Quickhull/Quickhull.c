#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define LEAD_RANK 0
#define STARTING_VECTOR_SIZE 64
#define CLOCK_RATE 512000000.0 // for clock to seconds conversion

typedef unsigned long long ticks;
static __inline__ ticks getticks(void) {
  unsigned int tbl, tbu0, tbu1;
  do {
       __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
       __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
       __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
     } while (tbu0 != tbu1);
  return ((((unsigned long long)tbu0) << 32) | tbl);
}

int myrank;
int numranks;

typedef struct Point {
    float x;
    float y;
} Point;

typedef struct Vector {
    Point * pts;
    size_t size;
    size_t cap;
} Vector;

typedef struct DistPoint {
    float dist;
    float x;
    float y;
} DistPoint;

typedef struct HullParams {
    Vector points;
    Point P;
    Point Q;
} HullParams;

void ensureReturnCode(int rc, char* point) {
    if (rc != MPI_SUCCESS) {
        char* message;
        int len;
        MPI_Error_string(rc, message, &len);
        printf("An error (%d) occured when %s\n", rc, point);
        printf("Message: %s\n", message);
        exit(-1);
    }
}

MPI_Datatype MPI_POINT;
void createPointType() {
    int blockLengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
    MPI_Aint offsets[2];

    offsets[0] = 0;
    offsets[1] = sizeof(float);

    MPI_Type_create_struct(2, blockLengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);
}

MPI_Datatype MPI_DIST_POINT;
void createDistPointType() {
    int blockLengths[3] = {1,1,1};
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Aint offsets[3];
    
    offsets[0] = 0;
    offsets[1] = sizeof(float);
    offsets[2] = 2*sizeof(float);

    MPI_Type_create_struct(3, blockLengths, offsets, types, &MPI_DIST_POINT);
    MPI_Type_commit(&MPI_DIST_POINT);
}

void readInData(char * fileName, Vector points, 
    int myrank, size_t stride) {
    MPI_File mpiFile;
    MPI_Status stat;

    int rc = MPI_File_open(MPI_COMM_WORLD, fileName, 
        MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
    ensureReturnCode(rc, "opening file");
    if (numPoints == 268435456) {
        rc = MPI_File_read_at(mpiFile, sizeof(Point) * myrank * stride, points.pts,
             (points.size / 2), MPI_POINT, &stat);
        ensureReturnCode(rc, "reading file");
        rc = MPI_File_read_at(mpiFile, 
            (myrank * stride * sizeof(Point)) + ((points.size / 2) * sizeof(Point)), 
            points + (points.size / 2),
            (points.size / 2), MPI_POINT, &stat);
    } else {
        rc = MPI_File_read_at(mpiFile, sizeof(Point) * myrank * stride, points.pts,
                  points.size, MPI_POINT, &stat);
    }
    ensureReturnCode(rc, "reading file");
    rc = MPI_File_close(&mpiFile);
    ensureReturnCode(rc, "closing file");
}

// returns the upper bound number from size given that
// the return value is a power of 2.
int sizeCap(size_t size) {
    size_t pref_size = pow(2, (int) (log(size)/log(2) + 1));
    size_t min_size = 16;
    return min_size < pref_size ? pref_size : min_size;
}

Vector vectorInit(size_t size) {
    Vector v;
    v.size = size;
    v.cap = sizeCap(size);
    v.pts = calloc(sizeCap(size), sizeof(Point));
    if (v.pts == NULL) {
        perror("COULD NOT ALLOC MEMORY SIZE");
    }
    return v;
}

void printVector(Vector v) {
    for (size_t i=0; i<v.size; i++) {
        printf("\tPoint: %f, %f\n", v.pts[i].x, v.pts[i].y);
    }
}

void freeVector(Vector v) {
    free(v.pts);
}

Point vectorPop(Vector * v) {
    (*v).size -= 1;
    Point p = (*v).pts[(*v).size];
    return p;
}

void vectorAdd(Vector * v, Point a) {
    size_t newSize = sizeCap((*v).cap + 1);
    (*v).size += 1;
    if ((*v).cap <= (*v).size) {
        Point * new_pts = realloc((*v).pts, newSize * sizeof(Point));
        if (new_pts != NULL) {
            (*v).pts = new_pts;
            (*v).cap = newSize;
        } else {
            perror("COULD NOT ALLOC MEMORY SIZE");
        }
    }
    (*v).pts[(*v).size - 1] = a;
}

// Vector concatHulls(Vector p1, Vector p2) {
//     Vector concat;
//     concat.size = p1.size + p2.size;
//     concat.pts = calloc(sizeCap(concat.size), sizeof(Point));

//     for (int p = 0; p < p1.size; ++p)
//         concat.pts[p] = p1.pts[p];
//     for (int p = 0; p < p2.size; ++p)
//         concat.pts[p + p1.size] = p2.pts[p];

//     freeVector(p1);
//     freeVector(p2);

//     return concat;
// }

float dist(Point a, Point b, Point c) {
    float tp = (c.x - b.x) * (b.y - a.y) - (b.x - a.x) * (c.y - b.y);
    float bt = pow(pow(c.x - b.x, 2.0) + pow(c.y - b.y, 2.0), 0.5);
    if (tp<0)
        return -tp/bt;
    return tp/bt;
}

bool same(Point a, Point b) {
    return a.x == b.x && a.y == b.y;
}

bool side(Point n, Point a, Point b) {
    return ((n.x-a.x) * (b.y-a.y) - (n.y-a.y) * (b.x-a.x)) < 0;
}

Point max(Vector v) {
    Point m = v.pts[0];
    for (size_t i = 1; i < v.size; i++) {
        Point cmp = v.pts[i];
        if (cmp.x > m.x || (cmp.x == m.x && cmp.y > m.y)) {
            m = cmp;
        }
    }
    return m;
}

Point min(Vector v) {
    Point m = v.pts[0];
    for (size_t i = 1; i < v.size; i++) {
        Point cmp = v.pts[i];
        if (cmp.x < m.x || (cmp.x == m.x && cmp.y < m.y)) {
            m = cmp;
        }
    }
    return m;
}

// void prepareToSendToLeadRank(Point * myPoints, int numPoints, int myrank, int numranks) {
//     // send and block on small size
//     MPI_Send(&numPoints, 1, MPI_UNSIGNED, LEAD_RANK, 0, MPI_COMM_WORLD);
//     MPI_Send(myPoints, numPoints, MPI_POINT, LEAD_RANK, 0, MPI_COMM_WORLD);
//     MPI_Send(&numPoints, 1, MPI_UNSIGNED, LEAD_RANK, 0, MPI_COMM_WORLD);
//     MPI_Send(myPoints, numPoints, MPI_POINT, LEAD_RANK, 0, MPI_COMM_WORLD);
// }

void mpi_max_func(void* inv, void* inoutv, int* len, MPI_Datatype* datatype) {
    Point* in = (Point*) inv;
    Point* inout = (Point*) inoutv;
    for (int i = 0; i < *len; i++) {
        if (in[i].x > inout[i].x || (in[i].x == inout[i].x && in[i].y > inout[i].y)) {
            inout[i] = in[i];
        }
    }
}

void mpi_dist_max_func(void* inv, void* inoutv, int* len, MPI_Datatype* datatype) {
    DistPoint* in = (DistPoint*) inv;
    DistPoint* inout = (DistPoint*) inoutv;
    for (int i = 0; i < *len; i++) {
        if (in[i].dist > inout[i].dist && in[i].dist > 0) {
            inout[i] = in[i];
        }
    }
}

void mpi_min_func(void* inv, void* inoutv, int* len, MPI_Datatype* datatype) {
    Point* in = (Point*) inv;
    Point* inout = (Point*) inoutv;
    for (int i = 0; i < *len; i++) {
        if (in[i].x < inout[i].x || (in[i].x == inout[i].x && in[i].y < inout[i].y)) {
            inout[i] = in[i];
        }
    }
}

void partition_set(Vector orig, Vector * a_v, Vector * b_v, Point a, Point b) {
    Point * pt = orig.pts;
    for (size_t i=0; i<orig.size; i++, pt++) {
        if (!same(*pt, a) && !same(*pt, b)) {
            if (side(*pt, a, b)) {
                vectorAdd(a_v, *pt);
            } else {
                vectorAdd(b_v, *pt);
            }
        }  
    }
}

// void filter_set(Vector * orig, Point a, Point b) {
//     for (size_t i=0; i<(*orig).size; i++) {
//         Point pt = (*orig).pts[i];
//         if (!same(pt, a) && !same(pt, b)) {
//             if (side(pt, a, b)) {
//                 vectorAdd(orig, pt);
//             }
//         }  
//     }
// }

Point reduceMin(Point m) {
    Point r;
    MPI_Op min_op;
    MPI_Op_create(mpi_min_func, 1, &min_op);
    MPI_Allreduce(&m, &r, 1, MPI_POINT, min_op, MPI_COMM_WORLD);
    return r;
}

Point reduceMax(Point m) {
    Point r;
    MPI_Op max_op;
    MPI_Op_create(mpi_max_func, 1, &max_op);
    MPI_Allreduce(&m, &r, 1, MPI_POINT, max_op, MPI_COMM_WORLD);
    return r;
}

void reduceDist(DistPoint * p, int len) {
    DistPoint * r = calloc(len, sizeof(DistPoint));
    MPI_Op op;
    MPI_Op_create(mpi_dist_max_func, 1, &op);
    MPI_Allreduce(p, r, len, MPI_DIST_POINT, op, MPI_COMM_WORLD);
    for (size_t i=0; i<len; i++) {
        p[i] = r[i];
    }
    free(r);
}

void outputHull(Vector hull) {
    for (int p = 0; p < hull.size; ++p) {
        printf("%.4f, %.4f\n", hull.pts[p].x, hull.pts[p].y);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        printf("FORMAT: %s <num points> <input file>\n", argv[0]);
        return false;
    }

    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);

    ticks start_ticks = getticks();
    ticks comm_overhead = 0;

    createPointType();
    createDistPointType();

    size_t totalPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &totalPoints))
        return false;

    if (myrank == 0) {
        printf("Number elements is %zu\n", totalPoints);
        printf("numranks is %d\n", numranks);
    }
    size_t numPoints = totalPoints / numranks;
    // add stragglers to the last rank
    if (numranks == myrank + 1) {
        numPoints += totalPoints % numPoints;
    }
    Vector points = vectorInit(numPoints);
    readInData(argv[2], points, myrank, totalPoints / numranks);

    Point l = min(points);
    Point r = max(points);
    ticks comm_start = getticks();
    l = reduceMin(l);
    r = reduceMax(r);
    comm_overhead += (getticks() - comm_start);


    Vector hull = vectorInit(2);
    hull.pts[0] = l;
    hull.pts[1] = r;

    #ifdef DEBUG
    if (myrank == 0) {
        printf("\tNew hull point: {%.2f, %.2f}\n", l.x, l.y);
        printf("\tNew hull point: {%.2f, %.2f}\n", r.x, r.y);
    }
    #endif

    Vector ls = vectorInit(0);
    Vector rs = vectorInit(0);

    partition_set(points, &ls, &rs, l, r);
    freeVector(points);
    HullParams * point_sets = calloc(2, sizeof(HullParams));
    point_sets[0].points = ls;
    point_sets[0].P = l;
    point_sets[0].Q = r;
    point_sets[1].points = rs;
    point_sets[1].P = r;
    point_sets[1].Q = l;

    int active = 1;
    int i = 0;
    size_t prev_size = 2;
    while (active) {
        #ifdef DEBUG
        if (myrank == 0) {
            printf("-- PASS %d --\n", i);
        }
        #endif
        HullParams * next_sets = calloc(prev_size * 2, sizeof(HullParams));
        DistPoint max_pts[prev_size];
        active = 0;
        for (int j=0; j<prev_size; j++) {
            #ifdef DEBUG
            if (myrank == 0) {
                printf("\tj=%d\n", j);
            }
            #endif
            Vector pts = point_sets[j].points;
            size_t num_pts = pts.size;
            if (num_pts == 0 || pts.pts == NULL) {
                HullParams h1;
                HullParams h2;
                h1.points.pts = NULL;
                h2.points.pts = NULL;
                next_sets[j*2] = h1;
                next_sets[j*2+1] = h2;
                max_pts[j].dist = -1;
                continue;
            }
            #ifdef DEBUG
            if (myrank == 0) {
                printf("\tj=%d: performing split\n", j);
            }
            #endif
            active = 1;
            Point P = point_sets[j].P;
            Point Q = point_sets[j].Q;
            
            DistPoint * dpts = calloc(num_pts, sizeof(DistPoint));
            dpts[0].dist = dist(pts.pts[0], P, Q);
            size_t max_ind = 0;
            for (int k=0; k<num_pts; k++) {
                dpts[k].dist = dist(pts.pts[k], P, Q);
                dpts[k].x = pts.pts[k].x;
                dpts[k].y = pts.pts[k].y;
                if (dpts[max_ind].dist < dpts[k].dist) {
                    max_ind = k;
                }
            }
            max_pts[j] = dpts[max_ind];
            free(dpts);
        }
        ticks comm_start = getticks();
        // SYNCHRONIZE MAX PT AND LOOP PARAMS
        // first, evaluate if we need to continue the recursion
        int act_result = 0;
        MPI_Allreduce(&active, &act_result, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        active = act_result;
        // perform reduction to find largest dist point
        reduceDist(max_pts, prev_size);
        comm_overhead += (getticks() - comm_start);
        
        for (int j=0; j<prev_size; j++) {
            Vector pts = point_sets[j].points;
            size_t num_pts = pts.size;

            Point P = point_sets[j].P;
            Point Q = point_sets[j].Q;
            Point c;
            c.x = max_pts[j].x;
            c.y = max_pts[j].y;
            #ifdef DEBUG
            if (myrank == 0 && max_pts[j].dist > 0) {
                printf("\tNew hull point: {%.2f, %.2f}\n", c.x, c.y);
                vectorAdd(&hull, c);
            }
            #endif
            if (num_pts == 0 || pts.pts == NULL) {
                continue;
            }
            
            Vector S0 = vectorInit(0);
            Vector S1 = vectorInit(0);
            Vector S2 = vectorInit(0);
            partition_set(pts, &S1, &S0, P, c);
            freeVector(pts);
            pts = vectorInit(0);
            partition_set(S0, &S2, &pts, c, Q);
            freeVector(S0);
            freeVector(pts);

            HullParams H1;
            HullParams H2;

            H1.points = S1;
            H2.points = S2;

            H1.P = P;
            H1.Q = c;
            H2.P = c;
            H2.Q = Q;

            next_sets[j*2] = H1;
            next_sets[j*2+1] = H2;
        }
        free(point_sets);
        point_sets = next_sets;
        prev_size *= 2;
        i++;
    }
    if (myrank == 0) {
        #ifdef DEBUG
        printf("\n\nTHE CONVEX HULL IS:\n");
        outputHull(hull);
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0) {
        ticks end_ticks = getticks();
        printf("Total Computation Time: %f\n", (double) (end_ticks - start_ticks) / CLOCK_RATE);
        printf("Communication Overhead: %f\n", (double) (comm_overhead) / CLOCK_RATE);
    }

    return true;
}
