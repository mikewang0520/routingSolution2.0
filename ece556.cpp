// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include <string>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cstdlib>
#include "ece556.h"

using namespace std;

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

  // allocate space for edgeCaps, edgeUtils, and edgeHistories
  rst->edgeCaps = (int*) malloc(rst->numEdges * sizeof(int));
  memset(rst->edgeCaps, rst->cap, rst->numEdges); // fill default capacities
  rst->edgeUtils = (int*) calloc(rst->numEdges, sizeof(int));
  rst->edgeHistories = (int*) calloc(rst->numEdges, sizeof(int));

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
    getline(linestream, item, ' '); // token 1 (p1y)
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
    getline(linestream, item, ' '); // token 4 (newCap)
    int newCap = stoi(item);
    rst->edgeCaps[edgeID] = newCap;

    // get ready to read next line
    linestream.str("");
    linestream.clear();
  }
  
  // clean up and return
  //myfile.close();
  return 1;
}

// This function creates a routing solution
int solveRouting(routingInst *rst)
{
  /*********** TO BE FILLED BY YOU **********/
  
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
	if (edgeID == -1) return -1;	
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
  return 1;
}


// Perform RRR on the given routing instance
int RRR(routingInst *rst, int useNetD, int useNetO) {
  int *netOrder = (int*) malloc(rst->numNets * sizeof(int));

  // determine net ordering
  if (useNetO) {
    netOrder = getNetOrder(rst);
  }
  else {
    for (int i = 0; i < rst->numNets; ++i) {
      // SHOULD THIS ONLY INCLUDE NETS OF NON-ZERO COST??
      netOrder[i] = rst->nets[i].id;
    }
  }

  // use net decomposition or not
  if (useNetD) {
    decomp(rst, netOrder);
  }
  else {
    // use algorithm from project part 1... how??
      // edit arguments to solveRouting to include netOrder?
  }

  // clean up and return
  free(netOrder);
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
