// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

#define DEBUG 1
#define ASCII_ZERO 48

// CITED RESOURCES:
// Command Line Parsing - https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html

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
    printf("Usage : ./ROUTE.exe <input_benchmark_name> <output_file_name> \n");
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
  
  // generate initial solution
  status = solveRouting(rst);
  if(status==0){
    printf("ERROR: running routing \n");
    release(rst);
    return 1;
  }
  
  // perform RRR
  if (useNetD !=0 || useNetO != 0) {
    status = -1;
    
    do {
      status = RRR(rst, useNetD, useNetO);
    } while (status != 0);
    
    if (status == -1) {
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
