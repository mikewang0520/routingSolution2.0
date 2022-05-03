// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.


#ifndef ECE556_H
#define ECE556_H

#include <stdio.h>
#include <cstdlib>
#include <cstring>
/**
 * A structure to represent a 2D Point. 
 */
typedef struct
{
  int x ; /* x coordinate ( >=0 in the routing grid)*/
  int y ; /* y coordinate ( >=0 in the routing grid)*/
  
} point ;


/**
 * A structure to represent a segment
 */
struct segment
{
  point p1 ; 	/* start point of a segment */
  point p2 ; 	/* end point of a segment */
  
  int numEdges ; 	/* number of edges in the segment*/
  int *edges ;  	/* array of edges representing the segment*/
  segment& operator = (const segment& other){
    this->p1 = other.p1;
    this->p2 = other.p2;
    this->numEdges = other.numEdges;
    free(this->edges);
    this->edges = (int*)malloc(other.numEdges*sizeof(int));
    memcpy(this->edges, other.edges, other.numEdges*sizeof(int));
    return *this;
  }
} ;


/**
 * A structure to represent a route
 */
typedef struct
{
  int numSegs ;  	/* number of segments in a route*/
  segment *segments ;  /* an array of segments (note, a segment may be flat, L-shaped or any other shape, based on your preference */
  
} route ;


/**
 * A structure to represent nets
 */
typedef struct
{
  
  int id ; 		/* ID of the net */
  int numPins ; 		/* number of pins (or terminals) of the net */
  point *pins ; 		/* array of pins (or terminals) of the net. */
  route nroute ;		/* stored route for the net. */

  int cost;
  
} net ;

/**
 * A structure to represent the routing instance
 */
typedef struct
{
  int gx ;		/* x dimension of the global routing grid */
  int gy ;		/* y dimension of the global routing grid */
  
  int cap ;
  
  int numNets ;	/* number of nets */
  net *nets ;		/* array of nets */
  
  int numEdges ; 	/* number of edges of the grid */
  int *edgeCaps; 	/* array of the actual edge capacities after considering for blockages */
  int *edgeUtils;	/* array of edge utilizations */
  int *edgeHistories;   /* index = edgeID, value = history value */
  int *edgeWeights;     /* index = edgeID, value = edge weight */
  
} routingInst ;


/* int readBenchmark(const char *fileName, routingInst *rst)
   Read in the benchmark file and initialize the routing instance.
   This function needs to populate all fields of the routingInst structure.
   input1: fileName: Name of the benchmark input file
   input2: pointer to the routing instance
   output: 1 if successful
*/
int readBenchmark(const char *fileName, routingInst *rst);


/* int solveRouting(routingInst *rst)
   This function creates a routing solution.
   
   input: pointer to the routing instance
   output: 1 if successful, 0 otherwise (e.g. the data structures are not populated) 
*/
int solveRouting(routingInst *rst);

/* in getEdgeWeight(routingInst *rst, int edgeID)
   Calculates the weight of an edge given its ID
   input1: pointer to the routing instance
   input2: edge ID
   output: weight of the edge
*/
int getEdgeWeight(routingInst *rst, int edgeID);


/* int getSegWeight(routingInst *rst, segment currSeg)
   This function calculates the cost of a segment by
   summing the weight of each edge.
   input1: pointer to the routing instance
   input2: current segment
   output: cost of the current segment
*/
int getSegWeight(routingInst *rst, segment &currSeg);


/* int getNetCost(routingInst *rst, net currNet)
   This function calculates the cost of a single net.
   input1: pointer to the routing instance
   input2: a net (NOT a pointer!!)
   output: cost of the net
*/
int getNetCost(routingInst *rst, net &currNet);


/* int getTotalCost(routingInst *rst)
   This function calculates the total cost of the
   routing solution.
   input1: pointer to the routing instance
   output: total cost of the routing instance
*/
int getTotalCost(routingInst *rst);


/* int* getNetOrder(routingInst *rst)
   This function calculates net costs that need to be
   ripped up and re-routed (i.e. cost != 0) and returns
   an ORDERED int array of the nets that will be re-routed.
   input1: pointer to the routing instance
   output: ORDERED int array of nets that will be re-routed
*/
int* getNetOrder(routingInst *rst);


/* void decomp(routingInst *rst)
   Net decomposition of all pins in all nets.
   input1: pointer to the routing instance
   output: doesn't return anything - rather, directy updates the nets and their pin ordering in rst
*/
void decomp(routingInst *rst);


/* int RRR(routingInst *rst, int useNetO, int iteration, int seed)
   Performs one iteration of "Rip-up and ReRoute".
   
   Calls "getNetOrder()" if using net ordering, otherwise
   uses the net order specified by the input file (possibly
   decomposed)
   
   input1: pointer to the routing instance
   input2: use net ordering if 1
   input3: current RRR iteration number
   input4: seed
   output: total cost of re-routed routing instance (use getTotalCost(rst))
*/
int RRR(routingInst *rst, int useNetO, int iteration, int seed);


/* int writeOutput(const char *outRouteFile, routingInst *rst)
   Write the routing solution obtained from solveRouting(). 
   Refer to the project link for the required output format.
   
   Finally, make sure your generated output file passes the evaluation script to make sure
   it is in the correct format and the nets have been correctly routed. The script also reports
   the total wirelength and overflow of your routing solution.
   
   input1: name of the output file
   input2: pointer to the routing instance
   output: 1 if successful, 0 otherwise 
*/
int writeOutput(const char *outRouteFile, routingInst *rst);


/* int release(routingInst *rst)
   Release the memory for all the allocated data structures. 
   Failure to release may cause memory problems after multiple runs of your program. 
   Need to recursively delete all memory allocations from bottom to top 
   (starting from segments then routes then individual fields within a net struct, 
   then nets, then the fields in a routing instance, and finally the routing instance)
   
   output: 1 if successful, 0 otherwise 
*/
int release(routingInst *rst);


#endif // ECE556_H
