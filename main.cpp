// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>

#define DEBUG 1
#define TIMESTAMPS 1
#define ASCII_ZERO 48
#define MAX_RUNTIME_MINS 15
#define MAX_PENALIZED_RUNTIME_MINS 20

// CITED RESOURCES:
// Command Line Parsing - https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html
// Elapsed Time - https://www.techiedelight.com/measure-elapsed-time-program-chrono-library/#:~:text=Since%20C%2B%2B11%2C%20the,It%20includes%20the%20%3Cchrono.

using namespace std;
int RRR_cycles;

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

  bool outOfTime = 0;
  bool stop_flag = 0;
  //if (useNetD != 0 || useNetO != 0) {
    printf("STARTING RRR...\n");
    auto start = chrono::steady_clock::now(); // starting time!
    auto max_runtime = chrono::steady_clock::now(); // max_runtime is based on the FIRST runtime. This might not always be true...

    RRR_cycles = 0;
    
    do {
      auto cycleStart = chrono::steady_clock::now(); // time of starting this cycle

      // print starting time
      if (TIMESTAMPS) printf("RRR Start Cycle %d: %ld minutes %ld seconds\n", RRR_cycles, chrono::duration_cast<chrono::seconds>(cycleStart - start).count() / 60, chrono::duration_cast<chrono::seconds>(cycleStart - start).count() % 60);

      // 15 minutes exceeded! try to run one more RRR then exit!
      if (chrono::duration_cast<chrono::seconds>(cycleStart - start).count() > (MAX_RUNTIME_MINS * 60)) {
	printf("RRR: 15 minutes exceeded! Attempting one more cycle...\n");
	outOfTime = 1;
      }

      // another cycle actually might exceed 20 minutes! Terminate now!
      if (RRR_cycles != 0) {
	if (chrono::duration_cast<chrono::seconds>(cycleStart - start).count() + chrono::duration_cast<chrono::seconds>(max_runtime - start).count() > ((MAX_PENALIZED_RUNTIME_MINS) * 60)) {
	  // another cycle might exceed 20 minutes! Terminate now!
	  printf("RRR: Another cycle might exceed 20 minutes! Terminating program...\n");
	  outOfTime = 1;
	  break;
	}
      }

      // ACTUALLY RUN RRR
      status = RRR(rst, useNetO, RRR_cycles);
      ++RRR_cycles;

      // get cycle end time
      auto cycleEnd = chrono::steady_clock::now();

      if (RRR_cycles == 1) {
	max_runtime = chrono::steady_clock::now();
      }

      if (TIMESTAMPS) {
	printf("RRR End: %ld minutes %ld seconds\n",
	       chrono::duration_cast<chrono::seconds>(cycleEnd - start).count() / 60,
	       chrono::duration_cast<chrono::seconds>(cycleEnd - start).count() % 60
	       );
      } 
      if ((useNetD == 0) || (useNetO == 0)){
        if (RRR_cycles > 4){
          stop_flag = 1;
        }
      }     
    } while (status > 0 && outOfTime == 0 && !stop_flag);
    
    if (status < 0) {
      printf("Num Cycles: %d\n", RRR_cycles);
      printf("ERROR: running RRR");
      release(rst);
      return 1;
    }
  //}

  if (outOfTime) printf("STOPPING RRR: Out Of Time!\n");
  
  // write the result
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
