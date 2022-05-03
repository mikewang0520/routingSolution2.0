// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include <string>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cstdlib>
#include <queue>
#include "ece556.h"

#define DEBUG 1      // general debug flag
#define EDGE_DEBUG 0 // used to turn edgeDataDump on or off
#define PROGRESS_DEBUG 1

using namespace std;

// FUNCTIONS

int getEdgeID(routingInst *rst, point p1, point p2)
{
  int x1;
  int y1;
  int x2;
  int y2;
  int gx = rst->gx;
  int gy = rst->gy;

  // assign point values such that p1 is "before" p2
  if (p2.x < p1.x || p2.y < p1.y)
  {
    // points were passed backwards
    x1 = p2.x;
    y1 = p2.y;
    x2 = p1.x;
    y2 = p1.y;
  }
  else
  {
    // points were passed properly
    x1 = p1.x;
    y1 = p1.y;
    x2 = p2.x;
    y2 = p2.y;
  }

  // ensure not the same point...
  if (x1 == x2 && y1 == y2)
  {
    cout << "getEdgeID: SAME POINT ERROR!" << endl;
    return -1;
  }
  // else if edge is horizontal
  else if (y1 == y2)
  {
    int edgeID = y2 * (gx - 1) + x2 -1;
    if (edgeID > rst->numEdges){
	return -1;
    }else{
      return edgeID;
    }
  }
  // else if edge is vertical
  else if (x1 == x2)
  {
    int edgeID = gx * gy + (y2 - 2) * gy +x2;
    if (edgeID > rst->numEdges){
	return -1;
    }else{
      return edgeID;
    }
  }
  // else not an edge
  else
  {
    cout << "getEdgeID: NOT AN EDGE!" << endl;
    return -1;
  }

  // default
  return -1;
}

void unpackEdgeID(routingInst *rst, int edgeID, point *p1, point *p2)
{
  // int x1 = p1.x;
  // int y1 = p1.y;
  // int x2 = p2.x;
  // int y2 = p2.y;
  int gx = rst->gx;
  int gy = rst->gy;
  int bias = ((gx - 1) * gy);

  // if edgeID vertical
  if (edgeID >= (gx - 1) * gy)
  {
    p2->y = (edgeID - bias) / gx + 1;
    p1->y = p2->y - 1;

    p2->x = (edgeID - bias) % gx;
    p1->x = p2->x;
  }
  // else edgeID horizontal
  else
  {
    p2->x = edgeID % (gx - 1) + 1;
    p1->x = p2->x - 1;

    p2->y = edgeID / (gx - 1);
    p1->y = p2->y;
  }

  return;
}

void dumpEdgeData(routingInst *rst)
{
  if (EDGE_DEBUG)
  {
    for (int i = 0; i < rst->numEdges; ++i)
    {
      printf("%d - cap: %d - util: %d - history: %d - weight: %d\n", i, rst->edgeCaps[i], rst->edgeUtils[i], rst->edgeHistories[i], rst->edgeWeights[i]);
    }
  }
}

int readBenchmark(const char *fileName, routingInst *rst)
{
  /*********** TO BE FILLED BY YOU **********/

  // open read file
  ifstream myfile;
  myfile.open(fileName, ios::in);

  // ensure file is open
  if (!myfile.is_open())
    return 0;

  // read x and y dimensions of global routing grid
  string line;
  getline(myfile, line);          // read line
  stringstream linestream(line);  // define linestream
  string item;                    // define
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
  getline(linestream, item, ' ');  // token 1 (cap)
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
  rst->nets = (net *)malloc(rst->numNets * sizeof(net));

  // read all pins of all nets (nested for loop)
  for (int i = 0; i < rst->numNets; ++i)
  {
    // iterate over all nets
    getline(myfile, line); // read line
    linestream.str(line);
    getline(linestream, item, ' ');   // token 0 (netName (e.g. "n0"))
    int netid = stoi(item.substr(1)); // extracts number from "n0"
    rst->nets[i].id = netid;
    getline(linestream, item, ' '); // token 1 (numPins)
    int numPins = stoi(item);
    rst->nets[i].numPins = numPins;
    linestream.str("");
    linestream.clear();

    // allocated per-net space for pins
    rst->nets[i].pins = (point *)malloc(rst->nets[i].numPins * sizeof(point));

    for (int j = 0; j < rst->nets[i].numPins; ++j)
    {
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
  rst->numEdges = ((rst->gy) * ((rst->gx) - 1)) + ((rst->gx) * ((rst->gy) - 1));

  // allocate space for edgeCaps, edgeUtils, edgeHistories, and edgeWeights
  rst->edgeCaps = (int *)malloc(rst->numEdges * sizeof(int)); // allocate capacities...
  for (int i = 0; i < rst->numEdges; ++i)
    rst->edgeCaps[i] = rst->cap;                                   // fill default capacities
  rst->edgeUtils = (int *)calloc(rst->numEdges, sizeof(int));      // allocate utilizations... (and set to 0)
  rst->edgeHistories = (int *)malloc(rst->numEdges * sizeof(int)); // allocate histories...
  for (int i = 0; i < rst->numEdges; ++i)
    rst->edgeHistories[i] = 1;                                  // default histories of 1 (makes sense AFTER "PART 1" in next step from main)
  rst->edgeWeights = (int *)calloc(rst->numEdges, sizeof(int)); // allocate edge weights... (and set to 0)

  if (DEBUG)
  {
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

  for (int i = 0; i < numBlockages; ++i)
  {
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

  if (DEBUG)
  {
    printf("After blockage reading...\n");
    dumpEdgeData(rst);
  }

  // clean up and return
  // myfile.close();
  return 1;
}

// This function creates a routing solution
int solveRouting(routingInst *rst)
{
  /*********** TO BE FILLED BY YOU **********/
  if (rst == 0)
    return -1; // failure!! (null rst)

  for (int i = 0; i < rst->numNets; ++i)
  {
    // iterate through all nets

    // there will be 1 less segment than there are pins
    // (e.g. 3 pins will only need 2 segments to connect)
    rst->nets[i].nroute.numSegs = rst->nets[i].numPins - 1;

    // allocate memory for segments
    rst->nets[i].nroute.segments = (segment *)malloc(rst->nets[i].nroute.numSegs * sizeof(segment));

    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
      // a segment is formed from a pin and its next pin in net
      point pin1 = rst->nets[i].pins[j];
      point pin2 = rst->nets[i].pins[j + 1];

      int xgap = abs(pin2.x - pin1.x);
      int ygap = abs(pin2.y - pin1.y);

      // allocate memory for minimum number of edges
      rst->nets[i].nroute.segments[j].numEdges = xgap + ygap;
      rst->nets[i].nroute.segments[j].edges = (int *)malloc(rst->nets[i].nroute.segments[j].numEdges * sizeof(int));

      int edgeCount = 0;

      // record start and end points of segments
      rst->nets[i].nroute.segments[j].p1 = pin1;
      rst->nets[i].nroute.segments[j].p2 = pin2;

      // pins have horizontal gap
      for (int k = 0; k < xgap; ++k)
      {
        // add horizontal edge to segment
        point currPoint;
        point nextPoint;

        currPoint.y = pin1.y;
        nextPoint.y = pin1.y;

        int edgeID = -1;

        // calculate edgeID
        if (pin2.x > pin1.x)
        {
          currPoint.x = pin1.x + k;
          nextPoint.x = pin1.x + k + 1;
          edgeID = getEdgeID(rst, currPoint, nextPoint);
        }
        else
        {
          currPoint.x = pin1.x - k;
          nextPoint.x = pin1.x - k - 1;
          edgeID = getEdgeID(rst, nextPoint, currPoint);
        }

        // store edge ID
        if (edgeID == -1)
          return -1; // failure!! (improper edge calculation)
        rst->nets[i].nroute.segments[j].edges[edgeCount] = edgeID;

        // increment edgeCount
        ++edgeCount;
      }

      // pins have vertical gap
      for (int k = 0; k < ygap; ++k)
      {
        // add vertical edge to segment
        point currPoint;
        point nextPoint;

        int edgeID = -1;

        currPoint.x = pin2.x;
        nextPoint.x = pin2.x;

        // calculate edgeID
        if (pin2.y > pin1.y)
        {
          currPoint.y = pin1.y + k;
          nextPoint.y = pin1.y + k + 1;
          edgeID = getEdgeID(rst, currPoint, nextPoint);
        }
        else
        {
          currPoint.y = pin1.y - k;
          nextPoint.y = pin1.y - k - 1;
          edgeID = getEdgeID(rst, nextPoint, currPoint);
        }

        // get and store edge ID
        if (edgeID == -1)
          return -1;
        rst->nets[i].nroute.segments[j].edges[edgeCount] = edgeID;

        // increment edgeCount
        ++edgeCount;
      }
    }
  }

  return 1; // success!
}

int getEdgeWeight(routingInst *rst, int edgeID)
{
  return rst->edgeWeights[edgeID];
}

int getSegWeight(routingInst *rst, segment &currSeg)
{
  int segWeight = 0;

  for (int i = 0; i < currSeg.numEdges; ++i)
  {
    segWeight += getEdgeWeight(rst, currSeg.edges[i]);
  }

  return segWeight;
}

int getNetCost(routingInst *rst, net &currNet)
{
  int netCost = 0;

  for (int i = 0; i < currNet.nroute.numSegs; ++i)
  {
    netCost += getSegWeight(rst, currNet.nroute.segments[i]);
  }

  return netCost;
}

int getTotalCost(routingInst *rst)
{
  int totalCost = 0;

  for (int i = 0; i < rst->numNets; ++i)
  {
    totalCost += getNetCost(rst, rst->nets[i]);
  }

  return totalCost;
}

void updateEdgeUtils(routingInst *rst)
{
  // set all edge utilizations to 0
  for (int i = 0; i < rst->numEdges; ++i)
    rst->edgeUtils[i] = 0;

  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {
    // for each segment
    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
      // for each edge
      for (int k = 0; k < rst->nets[i].nroute.segments[j].numEdges; ++k)
      {
        // increment edge's utilization
        rst->edgeUtils[rst->nets[i].nroute.segments[j].edges[k]] += 1;
      }
    }
  }
}

void updateEdgeWeights(routingInst *rst)
{
  // set all edge weights to 0
  for (int i = 0; i < rst->numEdges; ++i)
    rst->edgeWeights[i] = 0;

  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {
    // for each segment
    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
      // for each edge
      for (int k = 0; k < rst->nets[i].nroute.segments[j].numEdges; ++k)
      {
        // calculate and set each edge's weight in edgeWeights (yes there will likely be some repetition)

        int currEdgeID = rst->nets[i].nroute.segments[j].edges[k];
        int o_k1 = rst->edgeUtils[currEdgeID] - rst->edgeCaps[currEdgeID];
        if (o_k1 < 0)
          o_k1 = 0;

        // w_k = o_k1 * h_k
        rst->edgeWeights[currEdgeID] = o_k1 * rst->edgeHistories[currEdgeID];
      }
    }
  }
}

void updateEdgeHistories(routingInst *rst)
{
  // edgesUpdated is used to avoid multiple history updates
  int *edgesUpdated = (int *)calloc(rst->numEdges, sizeof(int)); // index = edgeID, value 1 = updated, value 0 = not updated

  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {
    // for each segment
    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
      // for each edge
      for (int k = 0; k < rst->nets[i].nroute.segments[j].numEdges; ++k)
      {
        // calculate and update each edge's history in edgeHistory (yes there will likely be some repetition)

        int currEdgeID = rst->nets[i].nroute.segments[j].edges[k];
        int o_k1 = rst->edgeUtils[currEdgeID] - rst->edgeCaps[currEdgeID];
        if (o_k1 < 0)
          o_k1 = 0;

        // increment history if overflow (AND IF EDGE HISTORY WAS NOT ALREADY UPDATED)
        if (edgesUpdated[currEdgeID] == 0 &&
            o_k1 > 0)
        {
          rst->edgeHistories[currEdgeID] += 1;
          edgesUpdated[currEdgeID] = 1; // "THIS HISTORY HAS NOW BEEN UPDATED"
        }
      }
    }
  }

  free(edgesUpdated);

  return;
}

int findDistance(point p1, point p2)
{
  return abs(p1.x - p2.x) + abs(p1.y - p2.y);
}

void findCorners(point p1, point p2, point *p3, point *p4)
{
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;

  // When p1 and p2 are on a straight line
  if (x2 == x1 || y2 == y1)
  {
    // p3 and p4 NULL means we have a straight line
    return;
  }
  // When p1 and p2 forms a rectangle
  else
  {
    p3->x = x1;
    p3->y = y2;
    p4->x = x2;
    p4->y = y1;
  }

  return;
}

void dumpNetPins(net currNet)
{
  printf("n%d\n", currNet.id);
  for (int i = 0; i < currNet.numPins; ++i)
  {
    printf("%d - (%d,%d)\n", i, currNet.pins[i].x, currNet.pins[i].y);
  }
}

// ...gross
void decomp(routingInst *rst)
{
  if (DEBUG)
  {
    printf("Start decomp...\n");
    for (int i = 0; i < rst->numNets; ++i)
    {
      dumpNetPins(rst->nets[i]);
    }
  }

  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {
    // allocate a new (decomposed) pins array
    // point *sortedPins = (point*) malloc(sizeof(rst->nets[i].pins));
    // memset(sortedPins, -1, rst->nets[i].numPins); // set all pin indexes to -1 by default

    // bubble sort (ish) the pins based on relevant distances (will run for index range [0,numPins-1] as it should)
    for (int j = 0; j < rst->nets[i].numPins; ++j)
    {
      int closestDistance = __INT_MAX__;
      int distance = __INT_MAX__;

      for (int k = j; k < rst->nets[i].numPins; ++k)
      {
        if (j == 0)
        {
          // find point closest to the origin to start
          point origin;
          origin.x = 0;
          origin.y = 0;

          if ((distance = findDistance(origin, rst->nets[i].pins[k])) < closestDistance)
          {
            std::swap(rst->nets[i].pins[0], rst->nets[i].pins[k]);
            closestDistance = distance;
          }
        }
        else if (j == 1)
        {
          // find closest pin to start pin (to form the first, ungodly rectangle...)
          if ((distance = findDistance(rst->nets[i].pins[0], rst->nets[i].pins[k])) < closestDistance)
          {
            std::swap(rst->nets[i].pins[1], rst->nets[i].pins[k]);
            closestDistance = distance;
          }
        }
        else
        {
          // ughhhh now we need to look at corners
          point p3;
          point p4;

          // initialization (to prevent "uninitialized" warnings)
          p3.x = -1;
          p3.y = -1;
          p4.x = -1;
          p4.y = -1;

          // propagate (maybe) p3 and p4 if corners are found. otherwise they'll both be NULL
          findCorners(rst->nets[i].pins[j - 2], rst->nets[i].pins[j - 1], &p3, &p4);

          // perform distance swaps from p1 and p2 (we don't know for sure if we have corners yet)

          // p1
          for (int l = j; l < rst->nets[i].numPins; ++l)
          {
            if ((distance = findDistance(rst->nets[i].pins[j - 2], rst->nets[i].pins[l])) < closestDistance)
            {
              std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
              closestDistance = distance;
            }
          }

          // p2
          for (int l = j; l < rst->nets[i].numPins; ++l)
          {
            if ((distance = findDistance(rst->nets[i].pins[j - 1], rst->nets[i].pins[l])) < closestDistance)
            {
              std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
              closestDistance = distance;
            }
          }

          // if false (p3.x is still -1), we can assume we have a straight line;
          // therefore, this function only performs distance swaps from p1 and p2,
          // so nothing more needs to happen and we're done for this loop cycle.
          if (p3.x != -1)
          {
            // perform distance swaps from p3 and p4 (rectangle corners)

            // p3
            for (int l = j; l < rst->nets[i].numPins; ++l)
            {
              if ((distance = findDistance(p3, rst->nets[i].pins[l])) < closestDistance)
              {
                std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
                closestDistance = distance;
              }
            }

            // p4
            for (int l = j; l < rst->nets[i].numPins; ++l)
            {
              if ((distance = findDistance(p4, rst->nets[i].pins[l])) < closestDistance)
              {
                std::swap(rst->nets[i].pins[j], rst->nets[i].pins[l]);
                closestDistance = distance;
              }
            }
          }
        }
      }
    }
  }

  if (DEBUG)
  {
    printf("End decomp...\n");
    for (int i = 0; i < rst->numNets; ++i)
    {
      dumpNetPins(rst->nets[i]);
    }
  }

  return;
}

void getNetOrder(routingInst *rst, int *netOrder)
{
  // propagate netOrder

  int netOrderIndex = 0;
  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {
    // if net has positive cost
    if (getNetCost(rst, rst->nets[i]) > 0)
    {
      // add net to netOrder
      netOrder[netOrderIndex++] = i;
    }
  }

  if (DEBUG)
  {
    printf("Pre getNetOrder Sorting...\n");
    printf("{");
    for (int i = 0; i < netOrderIndex; ++i)
    {
      printf("n%d = %d", netOrder[i], getNetCost(rst, rst->nets[netOrder[i]]));
      if (i + 1 != netOrderIndex)
        printf(", ");
    }
    printf("}\n");
  }

  int *netCostList = (int *)malloc(rst->numNets * sizeof(int));
  for (int i = 0; i < rst->numNets; ++i)
  {
    if (PROGRESS_DEBUG)
    {
      printf("Getting net cost of index %d...\r", i);
      fflush(stdout);
    }
    netCostList[i] = getNetCost(rst, rst->nets[i]);
  }
  if (PROGRESS_DEBUG)
  {
    printf("\n");
  }

  // reverse bubble sort netOrder (in descending order) -> O(n2)

  // netOrderIndex = num elements in netOrder
  for (int i = netOrderIndex - 1; i >= 0; --i)
  {
    if (PROGRESS_DEBUG)
    {
      printf("Now ordering index %d...\r", i);
      fflush(stdout);
    }
    for (int j = 0; j < i; ++j)
    {
      int jCost = netCostList[netOrder[j]];
      int j1Cost = netCostList[netOrder[j + 1]];
      if (jCost < j1Cost)
      {
        // swap elements (pushes higher cost towards left)
        std::swap(netOrder[j], netOrder[j + 1]);
      }
    }
  }
  if (PROGRESS_DEBUG)
    printf("\n");

  if (DEBUG)
  {
    printf("Post getNetOrder Sorting...\n");
    printf("{");
    for (int i = 0; i < netOrderIndex; ++i)
    {
      printf("n%d = %d", netOrder[i], getNetCost(rst, rst->nets[netOrder[i]]));
      if (i + 1 != netOrderIndex)
        printf(", ");
    }
    printf("}\n");
  }

  free(netCostList);
}

// Rips up all nets in the netOrder and updates edge utilizations therein
void RU(routingInst *rst, int netIndex)
{
  if (DEBUG)
  {
    printf("Before RU... \n");
    dumpEdgeData(rst);
  }


  // free all routing info for all nets in "netOrder" (until we hit a "net index" (NOT ID) of -1)
  // for each segment
  for (int j = 0; j < rst->nets[netIndex].nroute.numSegs; ++j)
  {
    // for each edge update its utilization
    for (int k = 0; k < rst->nets[netIndex].nroute.segments[j].numEdges; ++k)
    {
      rst->edgeUtils[rst->nets[netIndex].nroute.segments[j].edges[k]] -= 1; // decrement EDGE UTILIZATION
    }
    // free the segment's edges array
    free(rst->nets[netIndex].nroute.segments[j].edges);
    rst->nets[netIndex].nroute.segments[j].numEdges = 0; // all edges ripped up
  }

  // free the current net's route's segments array (NO NEED TO FREE SEGMENTS ARRAY -> NUMSEGS WILL ALWAYS BE THE SAME)
  // free(rst->nets[netOrder[i]].nroute.segments);

  // reset the current net's route's numSegs to 0 (NUMSEGS WILL ALWAYS BE THE SAME)
  // rst->nets[netOrder[i]].nroute.numSegs = 0;

  //updateEdgeUtils(rst); <- NOTE: you might need this
  if (DEBUG)
  {
    printf("After RU...\n");
    dumpEdgeData(rst);
  }

  return;
}

// void printQueue(queue<subpath> q)
// {
//   printf("Printing Queue...\n");
//   int qsize = (int)q.size();
//   for (int i = 0; i < qsize; i++)
//   {
//     subpath top = q.front();
//     printf("{%d,%d}{%d,%d} %d\n",
//            top.from.x,
//            top.from.y,
//            top.to.x,
//            top.to.y,
//            top.cost);
//     q.pop();
//   }
// }

// void reverseQueue(queue<subpath> q, queue<subpath> *rev_q)
// {
//   stack<subpath> st; // stack -> LIFO (last-in, first-out)

//   // put queue in stack
//   while (!q.empty())
//   {
//     st.push(q.front());
//     q.pop();
//   }

//   // put stack back in queue
//   while (!st.empty())
//   {
//     (*rev_q).push(st.top());
//     st.pop();
//   }
// }

int generateInitialSolution(routingInst *rst, int netIndex){
/*********** TO BE FILLED BY YOU **********/
  if (rst == 0)
    return -1; // failure!! (null rst)

  rst->nets[netIndex].nroute.numSegs = rst->nets[netIndex].numPins - 1;

  // allocate memory for segments
  rst->nets[netIndex].nroute.segments = (segment *)malloc(rst->nets[netIndex].nroute.numSegs * sizeof(segment));

  for (int j = 0; j < rst->nets[netIndex].nroute.numSegs; ++j)
  {
    // a segment is formed from a pin and its next pin in net
    point pin1 = rst->nets[netIndex].pins[j];
    point pin2 = rst->nets[netIndex].pins[j + 1];

    int xgap = abs(pin2.x - pin1.x);
    int ygap = abs(pin2.y - pin1.y);

    // allocate memory for minimum number of edges
    rst->nets[netIndex].nroute.segments[j].numEdges = xgap + ygap;
    rst->nets[netIndex].nroute.segments[j].edges = (int *)malloc(rst->nets[netIndex].nroute.segments[j].numEdges * sizeof(int));

    int edgeCount = 0;

    // record start and end points of segments
    rst->nets[netIndex].nroute.segments[j].p1 = pin1;
    rst->nets[netIndex].nroute.segments[j].p2 = pin2;

    // pins have horizontal gap
    for (int k = 0; k < xgap; ++k)
    {
      // add horizontal edge to segment
      point currPoint;
      point nextPoint;

      currPoint.y = pin1.y;
      nextPoint.y = pin1.y;

      int edgeID = -1;

      // calculate edgeID
      if (pin2.x > pin1.x)
      {
        currPoint.x = pin1.x + k;
        nextPoint.x = pin1.x + k + 1;
        edgeID = getEdgeID(rst, currPoint, nextPoint);
      }
      else
      {
        currPoint.x = pin1.x - k;
        nextPoint.x = pin1.x - k - 1;
        edgeID = getEdgeID(rst, nextPoint, currPoint);
      }

      // store edge ID
      if (edgeID == -1)
        return -1; // failure!! (improper edge calculation)
      rst->nets[netIndex].nroute.segments[j].edges[edgeCount] = edgeID;

      // increment edgeCount
      ++edgeCount;
    }

    // pins have vertical gap
    for (int k = 0; k < ygap; ++k)
    {
      // add vertical edge to segment
      point currPoint;
      point nextPoint;

      int edgeID = -1;

      currPoint.x = pin2.x;
      nextPoint.x = pin2.x;

      // calculate edgeID
      if (pin2.y > pin1.y)
      {
        currPoint.y = pin1.y + k;
        nextPoint.y = pin1.y + k + 1;
        edgeID = getEdgeID(rst, currPoint, nextPoint);
      }
      else
      {
        currPoint.y = pin1.y - k;
        nextPoint.y = pin1.y - k - 1;
        edgeID = getEdgeID(rst, nextPoint, currPoint);
      }

      // get and store edge ID
      if (edgeID == -1)
        return -1;
      rst->nets[netIndex].nroute.segments[j].edges[edgeCount] = edgeID;

      // increment edgeCount
      ++edgeCount;
    }
  }
  return 1; // success!
}

int RR(routingInst *rst, int netIndex, int iteration){
  srand(20);
	
  if(DEBUG){
    printf("Starting RR...\n");
  }

  int status = generateInitialSolution(rst, netIndex);
  if (status == -1){
    printf("Failed to generate initial solution in RR\n");
  }

  printf("Generated initial solution\n");


  int dx = iteration + 1;
  int dy = iteration + 1;
  point curr_point;
  point nxt_point;
  int edgeID;

  // iterate through the segments of the net -> store best solution in segment structure
  for(int i=0; i<rst->nets[netIndex].nroute.numSegs; ++i){
    point pin1 = rst->nets[netIndex].pins[i];
    point pin2 = rst->nets[netIndex].pins[i + 1];

    int xgap = abs(pin2.x - pin1.x);
    int ygap = abs(pin2.y - pin1.y);
    
    int numedges = xgap + ygap + (dx * 2);

    // NEED TO ALLOCATE MEMORY //
    iteration = iteration + 1;

    // initialize the temporary segment to compar to current solution
    segment new_seg;
    new_seg.p1 = pin1;
    new_seg.p2 = pin2;
    new_seg.edges = (int*)calloc(numedges, sizeof(int));

    // CASE 1 => horizontal line //
    if (ygap == 0){ 
      // (1.a) checking 'n' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((pin1.y + dy) <= (rst->gy - 1)){
        nxt_point = pin1;
        // add edges moving up by dy to new_seg
        for(int j=0; j<dy; ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move left or right adding edges until nxt_point.x == pin2.x
        if(pin2.x > pin1.x){ // if pin2 if right of pin1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is left of pin 1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving down by dy to new_seg
        for(int j=0; j<dy; ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
   
      // (1.b) checking 'u' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((pin1.y - dy) >= 0){
        nxt_point = pin1;
        // add edges moving down by dy to new_seg
        for(int j=0; j<dy; ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move left or right adding edges until nxt_point.x == pin2.x
        if(pin2.x > pin1.x){ // if pin2 if right of pin1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is left of pin 1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (getEdgeID(rst, curr_point, nxt_point) != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving up by dy to new_seg
        for(int j=0; j<dy; ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
    
    } 
    // CASE 2 => vertical line //
    else if (xgap == 0){ 
      // (2.a) checking 'c' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((pin1.x - dx) >= 0){
        nxt_point = pin1;
        // add edges moving left by dx to new_seg
        for(int j=0; j<dx; ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move up or down adding edges until nxt_point.y == pin2.y
        if(pin2.y > pin1.y){ // pin2 is above pin1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is below pin 1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx; ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
   
      // (2.b) checking horizontally mirrored 'c' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((pin1.x + dx) <= (rst->gx - 1)){
        nxt_point = pin1;
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx; ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move up or down adding edges until nxt_point.y == pin2.y
        if(pin2.y > pin1.y){ // pin2 is above pin1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is below pin 1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx; ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
    
    } 
    // CASE 3 => not a line //
    else { 
     // (3.a) checking 'c' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((min(pin1.x, pin2.x) - dx) >= 0){
        nxt_point = pin1;
        // add edges moving left by dx to new_seg
        for(int j=0; j<dx+max(0,(pin1.x-pin2.x)); ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move up or down adding edges until nxt_point.y == pin2.y
        if(pin2.y > pin1.y){ // pin2 is above pin1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is below pin 1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx+max(0,(pin2.x-pin1.x)); ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
   
 //     (3.b) checking horizontally mirrored 'c' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((max(pin1.x, pin2.x) + dx) <= (rst->gx - 1)){
        nxt_point = pin1;
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx+max(0,(pin2.x-pin1.x)); ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move up or down adding edges until nxt_point.y == pin2.y
        if(pin2.y > pin1.y){ // pin2 is above pin1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is below pin 1
          for(int j=0; j<ygap; ++j){
            curr_point = nxt_point;
            nxt_point.y = curr_point.y - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving right by dx to new_seg
        for(int j=0; j<dx+max(0,(pin1.x-pin2.x)); ++j){
          curr_point = nxt_point;
          nxt_point.x = curr_point.x - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
 
  //   (3.c) checking 'n' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((max(pin1.y, pin2.y) + dy) <= (rst->gy - 1)){
        nxt_point = pin1;
        // add edges moving up by dy to new_seg
        for(int j=0; j<dy+max(0,(pin2.y-pin1.y)); ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move left or right adding edges until nxt_point.x == pin2.x
        if(pin2.x > pin1.x){ // if pin2 if right of pin1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is left of pin 1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving down by dy to new_seg
        for(int j=0; j<dy+max(0,(pin1.y-pin2.y)); ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
           rst->nets[netIndex].nroute.segments[i] = new_seg;
        }
      }
    
      // (3.d) checking 'u' shaped segment //
      new_seg.numEdges = 0; // reset new_segment
      if ((min(pin1.y, pin2.y) - dy) >= 0){
        nxt_point = pin1;
        // add edges moving down by dy to new_seg
        for(int j=0; j<dy+max(0,(pin1.y-pin2.y)); ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y - 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // move left or right adding edges until nxt_point.x == pin2.x
        if(pin2.x > pin1.x){ // if pin2 if right of pin1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x + 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        } else { // pin2 is left of pin 1
          for(int j=0; j<xgap; ++j){
            curr_point = nxt_point;
            nxt_point.x = curr_point.x - 1;
            // check if valid edge
            edgeID = getEdgeID(rst, curr_point, nxt_point);
            if (edgeID != -1){
              // store edge and increment the total number edges in temp segment
              new_seg.edges[new_seg.numEdges] = edgeID;
              new_seg.numEdges = new_seg.numEdges + 1;
            } else {
              return -1; // tried to get invalid edge id --> out of bounds
            }
          }
        }
        // add edges moving up by dy to new_seg
        for(int j=0; j<dy+max(0,(pin2.y-pin1.y)); ++j){
          curr_point = nxt_point;
          nxt_point.y = curr_point.y + 1;
          // check if valid edge
          edgeID = getEdgeID(rst, curr_point, nxt_point);
          if (edgeID != -1){
            // store edge and increment the total number edges in temp segment
            new_seg.edges[new_seg.numEdges] = edgeID;
            new_seg.numEdges = new_seg.numEdges + 1;
          } else {
            return -1; // tried to get invalid edge id --> out of bounds
          }
        }
        // compare total weight of this segment to curr solution -> replace if better
        if(getSegWeight(rst, rst->nets[netIndex].nroute.segments[i]) > getSegWeight(rst, new_seg)){
          rst->nets[netIndex].nroute.segments[i] = new_seg;
          // iterate through edges in new_seg and increment ... ect
        }
      }
    }
  free(new_seg.edges);
  }
  return 0;
}

// Perform RRR on the given routing instance
int RRR(routingInst *rst, int useNetO, int iteration)
{
  // UPDATE EDGE HISTORIES HERE
  updateEdgeUtils(rst);     // "sets" on first iteration (calc u_k1)
  updateEdgeHistories(rst); // will always "update" (calc h_k)
  updateEdgeWeights(rst);   // "sets" on first iteration (calc w_k which is a function of u_k1 and h_k)

  // count how many non-zero nets there are
  int nonzeroNets = 0;
  for (int i = 0; i < rst->numNets; ++i)
  {
    if (getNetCost(rst, rst->nets[i]) > 0)
    {
      ++nonzeroNets;
    }
  }

  int netOrder[nonzeroNets];
  for (int i = 0; i < nonzeroNets; ++i)
    netOrder[i] = -1; // default all values to -1

  // determine net ordering (needs w_k !!)
  if (useNetO)
  {
    getNetOrder(rst, netOrder); // will only RRR nets with non-zero cost
  }
  else
  {
    int netOrderIndex = 0;
    for (int i = 0; i < rst->numNets; ++i)
    {
      if (getNetCost(rst, rst->nets[i]) > 0)
      {                                // NOTE: NOT THE MOST EFFICIENT
        netOrder[netOrderIndex++] = i; // will RRR all non-zero nets, in the order specified by the input file
      }
    }
  }

  // Rip Up
  for (int i = 0; i < nonzeroNets; ++i)
  {
    printf("RIPPING UP NET %d/%d... ID = %d\n", i, nonzeroNets-1, netOrder[i]);
    RU(rst, netOrder[i]);

    printf("REROUTING...\n");
    int check = RR(rst, netOrder[i], iteration);
    if (check == -1){
      printf("PROBLEM WITH RR -> RETURNED -1");
    }
    if (PROGRESS_DEBUG)
      printf("\n");

   // if (getTotalCost(rst) == 0)
   //   break;

    //updateEdgeUtils(rst); <- NOTE: this can be replaced
    //updateEdgeWeights(rst);
  }

  // UPDATE EDGE UTILS AFTER RU (weights should NOT be updated ??? otherwise w_k would be a function of o_k instead of o_k1
  // updateEdgeUtils(rst); // re-route will create new edges
  // updateEdgeWeights(rst);

  // Re Route
  // RR(rst, netOrder, nonzeroNets);

  // update values
  updateEdgeUtils(rst); //<- NOTE: keep track of while you optimize segments
  updateEdgeWeights(rst);

  int totalCost = getTotalCost(rst);
  printf("RR solution cost: %d\n", totalCost);

  // update edgeUtils and edgeWeights before returning a cost function (???)
  // updateEdgeUtils(rst);
  // updateEdgeWeights(rst);
  return totalCost;
}

// Write the routing solution
int writeOutput(const char *outRouteFile, routingInst *rst)
{
  // declare file stream variable
  ofstream fileOut;
  // open the output file
  fileOut.open(outRouteFile);

  // write net segments to fileOut
  for (int i = 0; i < rst->numNets; ++i) // enumerates through nets
  {
    fileOut << "n" << rst->nets[i].id << endl;

    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
      // iterates through endpoints of routes

      for (int k = 0; k < rst->nets[i].nroute.segments[j].numEdges; ++k)
      {
        point p1;
        point p2;

        unpackEdgeID(rst, rst->nets[i].nroute.segments[j].edges[k], &p1, &p2);

        fileOut << "(" << p1.x << "," << p1.y << ")-(" << p2.x << "," << p2.y << ")" << endl;
      }
    }

    fileOut << "!" << endl;
  }

  // close the output file
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
int release(routingInst *rst)
{
  /*********** TO BE FILLED BY YOU **********/
  if (!rst)
    return 0; // failure if rst is NULL

  // for each net
  for (int i = 0; i < rst->numNets; ++i)
  {

    // free all pins
    free(rst->nets[i].pins);

    // for each segment
    for (int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
    {
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
