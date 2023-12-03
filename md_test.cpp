/*
 * Program: 	md_test.cpp
 * Summary: 	This program is a proof of concept for a Molecular Dynamic model. This program
 *              only uses sequential programming techniques. This is a proof of concept.
 * Programmer:	Sean B. Higgins
 * Start Date:	November 27, 2023
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


// the matrix-width and the number of work-items per work-group:
// note: the matrices are actually NUM_ATOMS * NUM_ATOMS * NUM_ATOMS and the work group sizes are LOCALSIZExLOCALSIZE:

#ifndef NUM_ELEMS
#define NUM_ELEMS		32
#endif

#ifndef NUM_ATOMS
#define NUM_ATOMS		NUM_ELEMS * NUM_ELEMS * NUM_ELEMS
#endif

#ifndef NUM_CELLS
#define	NUM_CELLS		8000
#endif

// For use with the neighbor matrix W. Considering a perfectly regular crystal latice, each atom has at most
// 27 neighbors. This includes itself, and will need to be handled accordingly.
#ifndef NUM_NEIGHBORS
#define NUM_NEIGHBORS   27
#endif

//#define   CSV

// Multi-dimensional arrays and 1-D arrays to hold the data for the simulation.

int atomsPerCell[NUM_CELLS];    // Keeps track of how many atoms are in each cell (k_c for all c in [0, NUM_CELLS-1])
int cellOccupancyMatrix[NUM_CELLS][NUM_ATOMS];  // Cell-layer occupancy matrix H. Using NUM_ATOMS as N_l, the larget number
                                                // of atoms per layer. This is only for testing, and the matrix will need
                                                // to be allocated dynamically when the actual N_l can be found.
int neighborMatrix[NUM_NEIGHBORS];              // The neighbor matrix W. Keeps track of each atom's neighbors, including itself.
                                                // This is important for interaction calcuations (i.e., force and energy).


int main(int argc, char* argv) {
    // (N1): Assign atoms to cells, based on position: atom i is in cell c_i,
    //       with multiple atoms allowed per cell.


    // (N2): Assign atoms in each cell to layers, where layer l includes the l^th members
    //       of each cell: The i^th atom in a cell will be in layer l_i. Each cell c
    //       requires a cell occupancy counter k_c, and due to multiple occupancy,
    //       incrementing k_c_i requires an 'atomic' operation.


    // (N3): Determine the maximum number of layers required, N_l, by finding the maximum
    //       of all the k_c using a reduction algorithm (i.e., N_l = max{k_c} for c = [0, NUM_CELLS-1]).
    //       I'm considering a modified binary search algorithm.

    
    // (N4): Build the cell-layer occupancy matrix H by setting each cell H_(c_i,l_i) = i
    //       for each atom i; the row and colum indices of H_(c,l) specify the cell (c <= N_c)
    //       and layer (l <= N_l). REMEMBER: N_c = NUM_CELLS while N_l is found above in section (N3)
    //       This means that the matrix H will need to be dynamically allocated at run-time once we
    //       know what the maximum number of layers N_l will be.
    //       Alternatively, for testing, we can just allocate the array globaly for testing by considering
    //       the largest possible value for N_l= NUM_ATOMS.

    
    // (N5): Construct the neighbor matrix W; for each atom i there are two nested loops to access
    //       the neighbors i' = H_(c,l), first over the corresponding cells c of c_i, and then over
    //       the layers l <= k_c for that cell c.


    // Force and Energy Evaluation: Handled over Time Steps
    // (F1): For each atom i, accumulate the total force f_i and (optionally) interaction energy u_i
    //       by considering the subset of atoms i' = W_(m,i) for m < m_i that lie within interaction
    //       range.


    // (F2): Sum the individual u_i to obtain the total interaction energy U (technically 2U) using
    //       a reduction algorithm. Again, a modified binary algorithm would work.


    return 0;
}