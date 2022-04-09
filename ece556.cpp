// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include <string>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cstdlib>
#include "ece556.h"

using namespace std;

int getEdgeID(routingInst *rst, point p1, point p2) {
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;
  int gx = rst->gx;
  int gy = rst->gy;

  // if edge is horizontal
  if (y1 == y2) {
    return y2*(gx-1)+x2-1;
  }
  // else if edge is vertical
  else if (x1 == x2) {
    return gx*gy+(y2-2)*gy+x2;
  }
  // else not an edge
  else {
    cout << "NOT AN EDGE!" << endl;
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

  // BLOCKAGE READING GOES HERE (???)
  
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
	// add horizontal edge to segment
	 // add horizontal edge to segment
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
    /*
    char *string = (char*) malloc(1000 * sizeof(char));
    if (!string) {
      fputs ("ERR: memory allocation for string failed", stderr);
      exit (EXIT_FAILURE);
    }
    */

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

// Release memory used
// NOTE: do these need to be in for loops?
int release(routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/
  if (!rst) return 0; // failure if rst is NULL
  
  // for each net
  for (int i=0; i<rst->numNets; ++i) {

    // free all pins
    free(rst->nets[i].pins);

    // for each segment
    for (int j=0; j<rst->nets[i].nroute.numSegs; ++j) {
      // free edges
      free(rst->nets[i].nroute.segments[j].edges);
    }

    // free segments array in route
    free(rst->nets[i].nroute.segments);
  }

  // free nets
  free(rst->nets);

  // free edgeCaps and edgeUtils
  // TODO
  
  // release the routing instance itself
  free(rst);

  return 1; // success!
}
