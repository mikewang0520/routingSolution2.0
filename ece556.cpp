// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include <string>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cstdlib>
#include "ece556.h"
#include <queue>
#include <unordered_set>

#define DEBUG 1

using namespace std;

// STRUCTS

struct subpath{
    point from ;
    point to ;
    int cost ;
    bool operator<(const subpath &rhs) const
    {
        return cost > rhs.cost;
    }
};
      
// got this point sturcture online
struct Point
{
  int x;
  int y;

  Point() { }
  Point(int x, int y)
  {
    this->x = x;
    this->y = y;
  }
  
  bool operator==(const Point& otherPoint) const
  {
    if (this->x == otherPoint.x && this->y == otherPoint.y) return true;
    else return false;
  }

  struct HashFunction
  {
    size_t operator()(const Point& point) const
    {
      size_t xHash = std::hash<int>()(point.x);
      size_t yHash = std::hash<int>()(point.y) << 1;
      return xHash ^ yHash;
    }
  };
};

// FUNCTIONS

int getEdgeID(routingInst *rst, point p1, point p2) {
  int x1;
  int y1;
  int x2;
  int y2;
  int gx = rst->gx;
  int gy = rst->gy;

  // assign point values such that p1 is "before" p2
  if (p2.x < p1.x || p2.y < p1.y) {
    // points were passed backwards
    x1 = p2.x;
    y1 = p2.y;
    x2 = p1.x;
    y2 = p1.y;
  }
  else {
    // points were passed properly
    x1 = p1.x;
    y1 = p1.y;
    x2 = p2.x;
    y2 = p2.y;
  }
  
  // ensure not the same point...
  if (x1 == x2 && y1 == y2) {
    cout << "getEdgeID: SAME POINT ERROR!" << endl;
    return -1;
  }
  // else if edge is horizontal
  else if (y1 == y2) {
    return y2*(gx-1)+x2-1;
  }
  // else if edge is vertical
  else if (x1 == x2) {
    return gx*gy+(y2-2)*gy+x2;
  }
  // else not an edge
  else {
    cout << "getEdgeID: NOT AN EDGE!" << endl;
    return -1;
  }

  // default
  return -1;
}

void unpackEdgeID(routingInst *rst, int edgeID, point *p1, point *p2) {
  //int x1 = p1.x;
  //int y1 = p1.y;
  //int x2 = p2.x;
  //int y2 = p2.y;
  int gx = rst->gx;
  int gy = rst->gy;
  int bias = ((gx-1)*gy);
  
  // if edgeID vertical
  if (edgeID >= (gx-1)*gy) {
    p2->y = (edgeID-bias)/gx + 1;
    p1->y = p2->y - 1;

    p2->x = (edgeID-bias)%gx;
    p1->x = p2->x;
  }
  // else edgeID horizontal
  else {
    p2->x = edgeID%(gx-1) + 1;
    p1->x = p2->x - 1;

    p2->y = edgeID/(gx-1);
    p1->y = p2->y;
  }

  return;
}

void dumpEdgeData(routingInst *rst) {
  for (int i=0; i<rst->numEdges; ++i) {
    printf("%d - cap: %d - util: %d - history: %d - weight: %d\n",i,rst->edgeCaps[i],rst->edgeUtils[i],rst->edgeHistories[i], rst->edgeWeights[i]);
  }
}

int readBenchmark(const char *fileName, routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/

  // open read file
  ifstream myfile;
  myfile.open(fileName, ios::in);

  // ensure file is open
  if (!myfile.is_open()) return 0;

  // read x and y dimensions of global routing grid
  string line;
  getline(myfile, line); // read line
  stringstream linestream(line); // define linestream
  string item; // define
  getline(linestream, item, ' '); // token 0 ("grid")
  getline(linestream, item, ' '); // token 1 (gx)
  rst->gx = stoi(item);
  getline(linestream, item, ' '); // token 2 (gy)
  rst->gy = stoi(item);
  linestream.str("");
  linestream.clear();

  // read capacity number
  getline(myfile, line); // read line
  linestream.str(line);
  getline(linestream, item, '\t'); // token 0 ("capacity")
  getline(linestream, item, ' '); // token 1 (cap)
  rst->cap = stoi(item);
  linestream.str("");
  linestream.clear();

  // read numNets
  getline(myfile, line); // read line
  linestream.str(line);
  getline(linestream, item, ' '); // token 0 ("num")
  getline(linestream, item, ' '); // token 1 ("net")
  getline(linestream, item, ' '); // token 2 (numNets)
  rst->numNets = stoi(item);
  linestream.str("");
  linestream.clear();

  // allocate space for nets
  rst->nets = (net*) malloc(rst->numNets * sizeof(net));

  // read all pins of all nets (nested for loop)
  for (int i = 0; i < rst->numNets; ++i) {
    // iterate over all nets
    getline(myfile, line); // read line
    linestream.str(line);
    getline(linestream, item, ' '); // token 0 (netName (e.g. "n0"))
    int netid = stoi(item.substr(1)); // extracts number from "n0"
    rst->nets[i].id = netid;
    getline(linestream, item, ' '); // token 1 (numPins)
    int numPins = stoi(item);
    rst->nets[i].numPins = numPins;
    linestream.str("");
    linestream.clear();

    // allocated per-net space for pins
    rst->nets[i].pins = (point*) malloc(rst->nets[i].numPins * sizeof(point));
    
    for (int j = 0; j < rst->nets[i].numPins; ++j) {
      // iterate over all pins within a net
      getline(myfile, line); // read line
      linestream.str(line);
      getline(linestream, item, '\t'); // token 0 (x)
      rst->nets[i].pins[j].x = stoi(item);
      getline(linestream, item, '\n'); // token 1 (y)
      rst->nets[i].pins[j].y = stoi(item);
      linestream.str("");
      linestream.clear();
    }
  }

  // calculate and store numEdges
  rst->numEdges = ((rst->gy) * ((rst->gx)-1)) + ((rst->gx) * ((rst->gy)-1));

  // allocate space for edgeCaps, edgeUtils, edgeHistories, and edgeWeights
  rst->edgeCaps = (int*) malloc(rst->numEdges * sizeof(int));      // allocate capacities...
  for (int i=0; i<rst->numEdges; ++i) rst->edgeCaps[i] = rst->cap; // fill default capacities
  rst->edgeUtils = (int*) calloc(rst->numEdges, sizeof(int));      // allocate utilizations... (and set to 0)
  rst->edgeHistories = (int*) malloc(rst->numEdges * sizeof(int)); // allocate histories...
  for (int i=0; i<rst->numEdges; ++i) rst->edgeHistories[i] = 1;   // default histories of 1 (makes sense AFTER "PART 1" in next step from main)
  rst->edgeWeights = (int*) calloc(rst->numEdges, sizeof(int));    // allocate edge weights... (and set to 0)

  if (DEBUG) {
    printf("After struct init...\n");
    dumpEdgeData(rst);
  }
  
  // BLOCKAGE READING GOES HERE
  // numBlockages
  // p1x p1y p2x p2y newCap
  // p1x p1y p2x p2y newCap
  // ...

  // get number of blockages
  getline(myfile, line); // read line
  linestream.str(line);
  getline(linestream, item, ' '); // token 0 (numBlockages)
  int numBlockages = stoi(item);
  linestream.str("");
  linestream.clear();

  for (int i=0; i<numBlockages; ++i) {
    // define two points
    point p1;
    point p2;
    
    // read line
    getline(myfile, line); // read line
    linestream.str(line);
    
    // read p1
    getline(linestream, item, ' '); // token 0 (p1x)
    p1.x = stoi(item);
    getline(linestream, item, '\t'); // token 1 (p1y)
    p1.y = stoi(item);

    // read p2
    getline(linestream, item, ' '); // token 2 (p2x)
    p2.x = stoi(item);
    getline(linestream, item, ' '); // token 3 (p2y)
    p2.y = stoi(item);

    // determine which edge lies between p1 and p2
    int edgeID = -1;
    edgeID = getEdgeID(rst, p1, p2);
    
    // update capacity of determined edge
    getline(linestream, item, '\n'); // token 4 (newCap)
    int newCap = stoi(item);
    rst->edgeCaps[edgeID] = newCap;

    // get ready to read next line
    linestream.str("");
    linestream.clear();
  }

  if (DEBUG) {
    printf("After blockage reading...\n");
    dumpEdgeData(rst);
  }
  
  // clean up and return
  //myfile.close();
  return 1;
}

// This function creates a routing solution
int solveRouting(routingInst *rst)
{
  /*********** TO BE FILLED BY YOU **********/
  if (rst == 0) return -1; // failure!! (null rst)
  
  for(int i=0; i < rst->numNets; ++i){
    // iterate through all nets

    // there will be 1 less segment than there are pins
    // (e.g. 3 pins will only need 2 segments to connect)
    rst->nets[i].nroute.numSegs = rst->nets[i].numPins - 1;

    // allocate memory for segments
    rst->nets[i].nroute.segments = (segment*) malloc(rst->nets[i].nroute.numSegs * sizeof(segment));
    
    for(int j = 0; j < rst->nets[i].nroute.numSegs; ++j){
      // a segment is formed from a pin and its next pin in net
      point pin1 = rst->nets[i].pins[j];
      point pin2 = rst->nets[i].pins[j+1];
      
      int xgap = abs(pin2.x - pin1.x);
      int ygap = abs(pin2.y - pin1.y);
      
      // allocate memory for minimum number of edges
      rst->nets[i].nroute.segments[j].numEdges = xgap + ygap;
      rst->nets[i].nroute.segments[j].edges = (int*) malloc(rst->nets[i].nroute.segments[j].numEdges * sizeof(int));

      int edgeCount = 0;

      // record start and end points of segments
      rst->nets[i].nroute.segments[j].p1 = pin1;
      rst->nets[i].nroute.segments[j].p2 = pin2;
      
      // pins have horizontal gap
      for (int k = 0; k < xgap; ++k) {
	// add horizontal edge to segment
	point currPoint;
	point nextPoint;

	currPoint.y = pin1.y;
	nextPoint.y = pin1.y;

	int edgeID = -1;

	// calculate edgeID
	if (pin2.x > pin1.x) {
	  currPoint.x = pin1.x + k;
	  nextPoint.x = pin1.x + k + 1;
	  edgeID = getEdgeID(rst, currPoint, nextPoint);
	}
	else {
	  currPoint.x = pin1.x - k;
	  nextPoint.x = pin1.x - k - 1;
	  edgeID = getEdgeID(rst, nextPoint, currPoint);
	}

	// store edge ID
	if (edgeID == -1) return -1; // failure!! (improper edge calculation)
	rst->nets[i].nroute.segments[j].edges[edgeCount] = edgeID;
	
	// increment edgeCount
	++edgeCount;
      }
      
      // pins have vertical gap
      for (int k = 0; k < ygap; ++k) {
	// add vertical edge to segment
        point currPoint;
        point nextPoint;

	int edgeID = -1;
	
        currPoint.x = pin2.x;
        nextPoint.x = pin2.x;

	// calculate edgeID
        if (pin2.y > pin1.y) {
          currPoint.y = pin1.y + k;
          nextPoint.y = pin1.y + k + 1;
	  edgeID = getEdgeID(rst, currPoint, nextPoint);
        }
        else {
          currPoint.y = pin1.y - k;
          nextPoint.y = pin1.y - k - 1;
	  edgeID = getEdgeID(rst, nextPoint, currPoint);
        }

        // get and store edge ID
        if (edgeID == -1) return -1;
        rst->nets[i].nroute.segments[j].edges[edgeCount] = edgeID;

        // increment edgeCount
        ++edgeCount;
      }
    }
  }
  
  return 1; // success!
}

int getEdgeWeight(routingInst *rst, int edgeID) {
  return rst->edgeWeights[edgeID];
}

int getSegWeight(routingInst *rst, segment currSeg) {
  int segWeight = 0;

  for (int i=0; i<currSeg.numEdges; ++i) {
    segWeight += getEdgeWeight(rst, currSeg.edges[i]);
  }

  return segWeight;
}

int getNetCost(routingInst *rst, net currNet) {
  int netCost = 0;

  for (int i=0; i<currNet.nroute.numSegs; ++i) {
    netCost += getSegWeight(rst, currNet.nroute.segments[i]);
  }

  return netCost;
}

int getTotalCost(routingInst *rst) {
  int totalCost = 0;

  for (int i=0; i<rst->numNets; ++i) {
    totalCost += getNetCost(rst, rst->nets[i]);
  }

  return totalCost;
}

void updateEdgeUtils(routingInst *rst) {
  // set all edge utilizations to 0
  memset(rst->edgeUtils, 0, rst->numEdges);

  // for each net
  for (int i=0; i<rst->numNets; ++i) {
    // for each segment
    for (int j=0; j<rst->nets[i].nroute.numSegs; ++j) {
      // for each edge
      for (int k=0; k<rst->nets[i].nroute.segments[j].numEdges; ++k) {
	// increment edge's utilization
	rst->edgeUtils[rst->nets[i].nroute.segments[j].edges[k]] += 1;
      }
    }
  }
}

void updateEdgeWeights(routingInst *rst) {
  // set all edge weights to 0
  memset(rst->edgeWeights, 0, rst->numEdges);

  // for each net
  for (int i=0; i<rst->numNets; ++i) {
    // for each segment
    for (int j=0; j<rst->nets[i].nroute.numSegs; ++j) {
      // for each edge
      for (int k=0; k<rst->nets[i].nroute.segments[j].numEdges; ++k) {
	// calculate and set each edge's weight in edgeWeights (yes there will likely be some repetition)
	
	int currEdgeID = rst->nets[i].nroute.segments[j].edges[k];
	int o_k1 = rst->edgeUtils[currEdgeID] - rst->edgeCaps[currEdgeID];
	if (o_k1 < 0) o_k1 = 0;
	
	// w_k = o_k1 * h_k
	rst->edgeWeights[currEdgeID] = o_k1 * rst->edgeHistories[currEdgeID];
      }
    }
  }
}

void updateEdgeHistories(routingInst *rst) {
  // edgesUpdated is used to avoid 
  int *edgesUpdated = (int*) calloc(rst->numEdges, sizeof(int)); // index = edgeID, value 1 = updated, value 0 = not updated

  // for each net
  for (int i=0; i<rst->numNets; ++i) {
    // for each segment
    for (int j=0; j<rst->nets[i].nroute.numSegs; ++j) {
      // for each edge
      for (int k=0; k<rst->nets[i].nroute.segments[j].numEdges; ++k) {
        // calculate and update each edge's history in edgeHistory (yes there will likely be some repetition)
	
        int currEdgeID = rst->nets[i].nroute.segments[j].edges[k];
        int o_k1 = rst->edgeUtils[currEdgeID] - rst->edgeCaps[currEdgeID];
        if (o_k1 < 0) o_k1 = 0;
	
	// increment history if overflow (AND IF EDGE HISTORY WAS NOT ALREADY UPDATED)
	if (edgesUpdated[currEdgeID] == 0 &&
	    o_k1 > 0) {
	  rst->edgeHistories[currEdgeID] += 1;
	  edgesUpdated[currEdgeID] = 1; // "THIS HISTORY HAS NOW BEEN UPDATED"
	}
      }
    }
  }
  
  return;
}

int findDistance(point p1, point p2) {
  return abs(p1.x - p2.x) + abs(p1.y - p2.y);
}

void findCorners(point p1, point p2, point *p3, point *p4) {
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;

  // When p1 and p2 are on a straight line
  if (x2 == x1 || y2 == y1) {
    // p3 and p4 NULL means we have a straight line
    return;
  }
  // When p1 and p2 forms a rectangle
  else {
    p3->x = x1;
    p3->y = y2;
    p4->x = x2;
    p4->y = y1;
  }
  
  return;
}

void dumpNetPins(net currNet) {
  printf("n%d\n",currNet.id);
  for (int i=0; i<currNet.numPins; ++i) {
    printf("%d - (%d,%d)\n",i,currNet.pins[i].x, currNet.pins[i].y);
  }
}

// ...gross
void decomp(routingInst *rst) {
  if (DEBUG) {
    printf("Start decomp...\n");
    for (int i=0; i<rst->numNets; ++i) {
      dumpNetPins(rst->nets[i]);
    }
  }
  
  // for each net
  for (int i=0; i<rst->numNets; ++i) {
    // allocate a new (decomposed) pins array
    //point *sortedPins = (point*) malloc(sizeof(rst->nets[i].pins));
    //memset(sortedPins, -1, rst->nets[i].numPins); // set all pin indexes to -1 by default
    
    // bubble sort (ish) the pins based on relevant distances (will run for index range [0,numPins-1] as it should)
    for (int j=0; j<rst->nets[i].numPins; ++j) {
      int closestDistance = __INT_MAX__;
      int distance = __INT_MAX__;
      
      for (int k=j; k<rst->nets[i].numPins; ++k) {
	if (j == 0) {
	  // find point closest to the origin to start
	  point origin;
	  origin.x = 0;
	  origin.y = 0;
	  
	  if ((distance = findDistance(origin, rst->nets[i].pins[k])) < closestDistance) {
	    std::swap(rst->nets[i].pins[0], rst->nets[i].pins[k]);
	    closestDistance = distance;
	  }
	}
	else if (j == 1) {
	  // find closest pin to start pin (to form the first, ungodly rectangle...)
	  if ((distance = findDistance(rst->nets[i].pins[0], rst->nets[i].pins[k])) < closestDistance) {
	    std::swap(rst->nets[i].pins[1], rst->nets[i].pins[k]);
	    closestDistance = distance;
	  }
	}
	else {
	  // ughhhh now we need to look at corners
	  point p3;
	  point p4;

	  // initialization (to prevent "uninitialized" warnings)
	  p3.x = -1;
	  p3.y = -1;
	  p4.x = -1;
	  p4.y = -1;
	  
	  // propagate (maybe) p3 and p4 if corners are found. otherwise they'll both be NULL
	  findCorners(rst->nets[i].pins[j-2], rst->nets[i].pins[j-1], &p3, &p4);

	  // perform distance swaps from p1 and p2 (we don't know for sure if we have corners yet)

	  // p1
	  for (int l=j; l<rst->nets[i].numPins; ++l) {
	    if ((distance = findDistance(rst->nets[i].pins[j-2], rst->nets[i].pins[l])) < closestDistance) {
	      std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
	      closestDistance = distance;
	    }
	  }

	  // p2
	  for (int l=j; l<rst->nets[i].numPins; ++l) {
	    if ((distance = findDistance(rst->nets[i].pins[j-1], rst->nets[i].pins[l])) < closestDistance) {
	      std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
	      closestDistance = distance;
	    }
	  }

	  // if false (p3.x is still -1), we can assume we have a straight line;
	  // therefore, this function only performs distance swaps from p1 and p2,
	  // so nothing more needs to happen and we're done for this loop cycle.
	  if (p3.x != -1) {
	    // perform distance swaps from p3 and p4 (rectangle corners)

	    // p3
	    for (int l=j; l<rst->nets[i].numPins; ++l) {
	      if ((distance = findDistance(p3, rst->nets[i].pins[l])) < closestDistance) {
		std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
		closestDistance = distance;
	      }
	    }

	    // p4
	    for (int l=j; l<rst->nets[i].numPins; ++l) {
              if ((distance = findDistance(p4, rst->nets[i].pins[l])) < closestDistance) {
                std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
                closestDistance = distance;
              }
            }
	  }
	}
      }
    }
  }

  if (DEBUG) {
    printf("End decomp...\n");
    for (int i=0; i<rst->numNets; ++i) {
      dumpNetPins(rst->nets[i]);
    }
  }

  return;
}

void getNetOrder(routingInst *rst, int *netOrder) {
  // propagate netOrder

  int netOrderIndex = 0;
  // for each net
  for (int i=0; i<rst->numNets; ++i) {
    // if net has positive cost
    if (getNetCost(rst, rst->nets[i]) > 0) {
      // add net to netOrder
      netOrder[netOrderIndex++] = i;
    }
  }

  if (DEBUG) {
    printf("Pre getNetOrder Sorting...\n");
    printf("{");
    for (int i=0; i<netOrderIndex; ++i) {
      printf("n%d = %d",netOrder[i], getNetCost(rst, rst->nets[netOrder[i]]));
      if (i+1 != netOrderIndex) printf(", ");
    }
    printf("}\n");
  }

  // reverse bubble sort netOrder (in descending order) -> O(n2)

  // netOrderIndex = num elements in netOrder
  for (int i = netOrderIndex; i >= 0; --i) {
    for (int j = netOrderIndex - 1; j > netOrderIndex - i; --j) {
      int jCost = getNetCost(rst, rst->nets[netOrder[j]]);
      int j1Cost = getNetCost(rst, rst->nets[netOrder[j-1]]);
      if (jCost > j1Cost) {
	// swap elements (pushes higher cost towards left)
	std::swap(netOrder[j], netOrder[j-1]);
      }
    }
  }

  if (DEBUG) {
    printf("Post getNetOrder Sorting...\n");
    printf("{");
    for (int i=0; i<netOrderIndex; ++i) {
      printf("n%d = %d",netOrder[i], getNetCost(rst, rst->nets[netOrder[i]]));
      if (i+1 != netOrderIndex) printf(", ");
    }
    printf("}\n");
  }
}

// Rips up all nets in the netOrder and updates edge utilizations therein
void RU(routingInst *rst, int *netOrder, int n) {
  if (DEBUG) {
    printf("Before RU... \n");
    dumpEdgeData(rst);
  }
  
  // free all routing info for all nets in "netOrder" (until we hit a "net index" (NOT ID) of -1)
  for (int i=0; i<n; ++i) {
    // for each segment
    for (int j=0; j<rst->nets[netOrder[i]].nroute.numSegs; ++j) {
      // for each edge update its utilization
      for (int k=0; k<rst->nets[netOrder[i]].nroute.segments[j].numEdges; ++k) {
        rst->edgeUtils[rst->nets[netOrder[i]].nroute.segments[j].edges[k]] -= 1; // decrement EDGE UTILIZATION
      }
      // free the segment's edges array
      free(rst->nets[netOrder[i]].nroute.segments[j].edges);
    }

    // free the current net's route's segments array (removes numEdges too)
    free(rst->nets[netOrder[i]].nroute.segments);

    // reset the current net's route's numSegs to 0
    rst->nets[netOrder[i]].nroute.numSegs = 0;
  }

  if (DEBUG) {
    printf("After RU...\n");
    dumpEdgeData(rst);
  }
  
  return;
}

void printQueue(queue<subpath> q) {
  printf("Printing Queue...\n");
  int qsize = (int)q.size();
  for(int i = 0; i < qsize; i++) {
	subpath top = q.front();
    printf("{%d,%d}{%d,%d} %d\n", 
	       top.from.x,
		   top.from.y,
		   top.to.x,
		   top.to.y,
		   top.cost
		   );
	q.pop();	
  }
}

queue<subpath> reverseQueue(queue<subpath> q) {
  stack<subpath>st; // stack -> LIFO (last-in, first-out)
  
  // put queue in stack
  while(!q.empty()){
    st.push(q.front());
    q.pop();
  }
  
  // put stack back in queue
  while(!st.empty()){
    q.push(st.top());
    st.pop();
  }

  return q;
}

queue<point> reverseQueue(queue<point> q) {
  stack<point>st; // stack -> LIFO (last-in, first-out)
  
  // put queue in stack
  while(!q.empty()){
    st.push(q.front());
    q.pop();
  }
  
  // put stack back in queue
  while(!st.empty()){
    q.push(st.top());
    
    st.pop();
  }
  
  return q;
}

void RR(routingInst *rst, int *netOrder, int n) {
  // Luke & Mike do maze routing
  if (DEBUG) {
    printf("Starting RR...\n");
  }

  // for each net that needs to be rerouted
  for (int i=0; i<n; ++i) {
    // propagate edgeAddCost array
    int edgeAddCost[rst->numEdges];
    for (int j=0; j<rst->numEdges; ++j) {
      // if 0 cap blockage, NEVER consider this edge!!
      if (rst->edgeCaps[j] == 0) {
	    edgeAddCost[j] = 999999;
      }
      else {
	    // cost to traverse this edge = 1 + its weight
	    edgeAddCost[j] = 1 + rst->edgeWeights[j];
      }
    }

    if (DEBUG) {
      printf("edgeAddCosts...\n");
      for (int j=0; j<rst->numEdges; ++j) {
	    printf("%d - %d\n", j, edgeAddCost[j]);
      }
    }

    // solve routing (from each pin to the next one IN ORDER)
    for (int j=0; j<rst->nets[netOrder[i]].numPins - 1; ++j) {
      // obtain the pins to be routed together
      point source = rst->nets[netOrder[i]].pins[j];
      point target = rst->nets[netOrder[i]].pins[j+1];
      
      // use Dijkstra's algorithm to solve between the two points ???
      printf("Maze routing between {%d, %d} and {%d, %d}...\n", source.x, source.y, target.x, target.y);

      // nodes that have been visited where min-cost path from starting point is known
      std::queue<subpath> RP; // revisit path

      // nodes that have been visited but min-cost path from starting node has not been found
      std::priority_queue<subpath> PQ;

      // set of all points with 
      std::unordered_set<Point, Point::HashFunction> visited_points;
        
      // initialize point
      point curr_pnt = source;
      int curr_pnt_cost = 0;
      RP.push(subpath{source,source,curr_pnt_cost}); // add the starting point to the min-cost path queue
      visited_points.insert(Point{curr_pnt.x, curr_pnt.y});
      
	  // cost function for A*
	  int fn;

      while((curr_pnt.x != target.x) || (curr_pnt.y != target.y)){ //  if current point is target point -> exit and retrace min-path
        if ((curr_pnt.x != 0) && (visited_points.find(Point{curr_pnt.x - 1, curr_pnt.y}) == visited_points.end())){ // if not bottom edge & not a visitied node
          // add line left
          point px = {curr_pnt.x - 1, curr_pnt.y};
          //fn = findDistance(px, target) + (getEdgeWeight(rst, getEdgeID(rst, curr_pnt, px))+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
          fn = findDistance(px, target) + (edgeAddCost[getEdgeID(rst, curr_pnt, px)]+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
		  PQ.push(subpath{curr_pnt, px, fn});
        }
        if ((curr_pnt.x != rst->gx - 1) && (visited_points.find(Point{curr_pnt.x + 1, curr_pnt.y}) == visited_points.end())){ // if not top edge & not a visitied node
          // add line right
          point px = {curr_pnt.x + 1, curr_pnt.y};
          //fn = findDistance(px, target) + (getEdgeWeight(rst, getEdgeID(rst, curr_pnt, px))+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
          fn = findDistance(px, target) + (edgeAddCost[getEdgeID(rst, curr_pnt, px)]+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
		  PQ.push(subpath{curr_pnt, px, fn});
        }
        if ((curr_pnt.y != 0) && (visited_points.find(Point{curr_pnt.x, curr_pnt.y - 1}) == visited_points.end())){ // if not left edge & not a visitied node
          // add line below
          point py = {curr_pnt.x, curr_pnt.y - 1};
          //fn = findDistance(py, target) + (getEdgeWeight(rst, getEdgeID(rst, curr_pnt, py))+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
          fn = findDistance(py, target) + (edgeAddCost[getEdgeID(rst, curr_pnt, py)]+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
		  PQ.push(subpath{curr_pnt, py, fn});  
        }
        if ((curr_pnt.y != rst->gy - 1) && (visited_points.find(Point{curr_pnt.x, curr_pnt.y + 1}) == visited_points.end())){ // if not right edge & not a visitied node
          // add line above
          point py = {curr_pnt.x, curr_pnt.y + 1};
          //fn = findDistance(py, target) + (getEdgeWeight (rst, getEdgeID(rst, curr_pnt, py))+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
		  fn = findDistance(py, target) + (edgeAddCost[getEdgeID(rst, curr_pnt, py)]+curr_pnt_cost); // f(n) = h(n) + g(n)  -> A* cost function
          PQ.push(subpath{curr_pnt, py, fn}); 
        }

        // get min_edge off pq for next iteration
        subpath min_edge = PQ.top();
		PQ.pop();
        // add to revisit path queue
        RP.push(min_edge);
        // set new current point
        curr_pnt = min_edge.to;
        // reset curr_pnt_cost
        curr_pnt_cost = min_edge.cost;
        // add point with known min-cost to point set
        visited_points.insert(Point{curr_pnt.x, curr_pnt.y});
      }
      printQueue(RP);
    } // END OF INNER FOR LOOP //

    point source_point = rst->nets[netOrder[i]].pins[0];
    queue<point>path;
    queue<subpath> RP_comparator;
    
    RP = reverseQueue(RP);
    RP_comparator = RP;
    RP_comparator.pop(); // get it one iteration ahead
    
    // add first two elements
    path.push(RP.front().to);
    path.push(RP.front().from);
    
    int cntr = 1;
    while((RP_comparator.front().from.x != source_point.x) || (RP_comparator.front().from.y != source_point.y)){
      if((RP.front().from.x == RP_comparator.front().to.x) & (RP.front().from.y == RP_comparator.front().to.y)){
        while(cntr != 0){
          RP.pop();
          --cntr;
        }
        cntr = 1;
        RP_comparator.pop();
        path.push(RP.front().from);
      } else{
        RP_comparator.pop();
        ++cntr;
      } 
    }
    path.push(source_point);
    path = reverseQueue(path); // final path of net => yay!


    // NEED TO ADD NEW PATH TO THE DATA STRUCT //


  } // END OF OUTER FOR LOOP //

  return;
}

// Perform RRR on the given routing instance
int RRR(routingInst *rst, int useNetO) {
  // UPDATE EDGE HISTORIES HERE
  updateEdgeUtils(rst); // "sets" on first iteration (calc u_k1)
  updateEdgeHistories(rst); // will always "update" (calc h_k)
  updateEdgeWeights(rst); // "sets" on first iteration (calc w_k which is a function of u_k1 and h_k)

  // count how many non-zero nets there are
  int nonzeroNets = 0;
  for (int i=0; i<rst->numNets; ++i) {
    if (getNetCost(rst, rst->nets[i]) > 0) {
      ++nonzeroNets;
    }
  }

  int netOrder[nonzeroNets];
  for (int i=0; i<nonzeroNets; ++i) netOrder[i] = -1; // default all values to -1
  
  // determine net ordering (needs w_k !!)
  if (useNetO) {
    getNetOrder(rst, netOrder); // will only RRR nets with non-zero cost
  }
  else {
	int netOrderIndex = 0;
    for (int i = 0; i < rst->numNets; ++i) {
      if (getNetCost(rst, rst->nets[i]) > 0) {
	    netOrder[netOrderIndex++] = i; // will RRR all non-zero nets, in the order specified by the input file
      }
    }
  }

  printf("nonzeroNets = %d\n",nonzeroNets);
  for (int i=0; i<nonzeroNets; ++i) {
	printf("%d\n",netOrder[i]);
  }

  // Rip Up
  RU(rst, netOrder, nonzeroNets);

  // UPDATE EDGE UTILS AFTER RU (weights should NOT be updated ??? otherwise w_k would be a function of o_k instead of o_k1
  updateEdgeUtils(rst); // re-route will create new edges
  //updateEdgeWeights(rst);

  // Re Route
  RR(rst, netOrder, nonzeroNets);
  
  // clean up and return
  //free(netOrder);

  if (DEBUG) {
    printf("RR solution cost: %d\n", getTotalCost(rst));
  }
  
  return getTotalCost(rst);
}

// Write the routing solution
int writeOutput(const char *outRouteFile, routingInst *rst)
{
  //declare file stream variable
  ofstream fileOut;
  //open the output file
  fileOut.open(outRouteFile);


  // write net segments to fileOut
  for (int i = 0; i < rst->numNets; ++i) // enumerates through nets
  {
    fileOut << "n" << rst->nets[i].id << endl;
    
    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j) {
      // iterates through endpoints of routes
      
      for (int k = 0; k < rst->nets[i].nroute.segments[j].numEdges; ++k) {
	point p1;
	point p2;
	
	unpackEdgeID(rst, rst->nets[i].nroute.segments[j].edges[k], &p1, &p2);
	
	fileOut << "(" <<
	   p1.x << "," <<
	   p1.y << ")-(" <<
	   p2.x << "," <<
	   p2.y << ")" << endl;
      }
    }
    
    fileOut << "!" << endl;
  }
  
  //close the output file
  fileOut.close();
  return 1;
}

/** 
 * Need to recursively delete all memory allocations from 
 * bottom to top (e.g., starting from segments, then routes
 * then individual fields within a net struct, then the
 * nets, then the fields in a routing instance, and finally
 * the routing instance) 
 */
int release(routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/
  if (!rst) return 0; // failure if rst is NULL
  
  // for each net
  for (int i=0; i<rst->numNets; ++i) {

    // free all pins
    free(rst->nets[i].pins);

    // for each segment
    for (int j=0; j<rst->nets[i].nroute.numSegs; ++j) {
      // free each edge
      free(rst->nets[i].nroute.segments[j].edges);
    }

    // free segments array in route
    free(rst->nets[i].nroute.segments);
  }

  // free nets
  free(rst->nets);

  // free edgeCaps and edgeUtils
  free(rst->edgeCaps);
  free(rst->edgeUtils);

  // free edgeHistories
  free(rst->edgeHistories);
  
  // free the routing instance itself
  free(rst);

  return 1; // success!
}
