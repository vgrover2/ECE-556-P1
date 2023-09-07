// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <cstring>

int main(int argc, char **argv)
{

 	if(argc!=5){
 		printf("Usage : ./ROUTE.exe -d=0 -n=0 <input_benchmark_name> <output_file_name> \n");
 		return 1;
 	}

 	int status;
	int d = argv[1][strlen(argv[1]) - 1] - '0';
	if(d != 0 && d != 1) {
		printf("d must equal 0 or 1\n");
	}

 	int n = argv[2][strlen(argv[2]) - 1] - '0';
	if(n != 0 && n!= 1) {
		printf("d must equal 0 or 1\n");
	}

	char *inputFileName = argv[3];
 	char *outputFileName = argv[4];




 	/// create a new routing instance
 	routingInst *rst = new routingInst;
	
 	/// read benchmark
    
 	status = readBenchmark(inputFileName, rst);
 	if(status==0){
 		printf("ERROR: reading input file \n");
 		return 1;
 	}

	if(d == 1) {
		status = netDecomposition(rst);
	}
    
 	/// run actual routing
    
 	status = solveRouting(rst);
    
	if(n == 1) {
    	status = RRR(rst);
	}
 	if(status==0){
 		printf("ERROR: running routing \n");
 		release(rst);
 		return 1;
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


