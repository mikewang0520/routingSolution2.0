// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>

#define DEBUG 0
#define ASCII_ZERO 48

// CITED RESOURCES:
// Command Line Parsing - https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html
// Elapsed Time - https://www.techiedelight.com/measure-elapsed-time-program-chrono-library/#:~:text=Since%20C%2B%2B11%2C%20the,It%20includes%20the%20%3Cchrono.

using namespace std;

int
main(int argc, char **argv)
{
  int useNetD = 0;
  int useNetO = 0;
  char c;
  
  while ((c = getopt(argc, argv, "d:n:")) != -1) {
    switch (c)
      {
      case 'd':
	useNetD = (int) optarg[1] - ASCII_ZERO;
	break;
      case 'n':
	useNetO = (int) optarg[1] - ASCII_ZERO;
	break;
      case '?':
	if (optopt == 'd' || optopt == 'n')
	  printf("Option -%c requires an argument.\n", optopt);
	else if (isprint(optopt))
	  printf("Unknown option '-%c'.\n",optopt);
	else
	  printf("Unknown option character '\\x%x'.\n", optopt);
	return 1;
      default:
	abort();
      }
  }

  if (DEBUG) {
    printf("useNetD = %d\n", useNetD);
    printf("useNetO = %d\n", useNetO);
  }
  
  if(argc < 3){
    printf("Usage : ./ROUTE.exe <optional arguments> <input_benchmark_name> <output_file_name> \n");
    return 1;
  }
  
  int status;
  char *inputFileName = argv[argc-2];
  char *outputFileName = argv[argc-1];
  
  // create a new routing instance
  routingInst *rst = new routingInst;
  
  // read benchmark
  status = readBenchmark(inputFileName, rst);
  if(status==0){
    printf("ERROR: reading input file \n");
    return 1;
  }

  if (useNetD == 1) {
    decomp(rst);
  }
  
  // generate initial solution
  status = solveRouting(rst);
  if(status==0){
    printf("ERROR: running routing \n");
    release(rst);
    return 1;
  }
  
  // perform RRR

  // status codes
  // 0 = perfect solution
  // # = total cost of routing instance
  // run while not perfect solution and also not over 15 minutes

  // NOTE: ADD TIMER OUT HERE IN MAIN.CPP
  
  if (useNetD != 0 || useNetO != 0) {
    auto start = chrono::steady_clock::now(); // starting time!
    int over15mins = 0;
    
    do {
      auto now = chrono::steady_clock::now();
      if (chrono::duration_cast<chrono::seconds>(now - start).count() > 900)
	over15mins = 1; // 15 minutes exceeded! run one more loop then exit
	
      status = RRR(rst, useNetO);
    } while (status > 0 && over15mins == 0);
    
    if (status < 0) {
      printf("ERROR: running RRR");
      release(rst);
      return 1;
    }
  }
  
    /// write the result
  status = writeOutput(outputFileName, rst);
  if(status==0){
    printf("ERROR: writing the result \n");
    release(rst);
    return 1;
  }
  
  release(rst);
  printf("\nDONE!\n");	
  return 0;
}
