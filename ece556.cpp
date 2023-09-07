// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <math.h>

using namespace std;
using namespace std::chrono;

typedef pair<int, Node*> costpair;

high_resolution_clock::time_point start;

int getEdge(int x, int y, int xDim) {
    return y * (xDim - 1) + x;
}


int getEdgeBetween(int x1, int y1, int x2, int y2, int xDim, int yDim) {
    if(x1 == x2) {
    // vertical edge
        if(y1 < y2) {
            // total horiz edges + getEdge(x1, y1)
            return yDim * (xDim - 1) + getEdge(x1, y1, xDim) + y1;
        }
        else {
            return yDim * (xDim - 1) + getEdge(x2, y2, xDim) + y2;
        }
    }
    else {
    // horizontal edge
        if(x1 < x2) {
            return getEdge(x1, y1, xDim);
        }
        else {
            return getEdge(x2, y2, xDim);
        }
  }
}



/**
 * Helper function to extract int values of x and y from a string
 *
 * toReturn determines which word in myString to return
 */
int findXY(string myString, int toReturn) {
  int i;
  int numSpaces = 0;
  int firstEnd = 0;
  int secondStart = 0;
  int secondEnd = 0;
  int thirdStart = 0;
  int thirdEnd = 0;
  int fourthStart = 0;
  int fourthEnd = 0;
  int fifthStart = 0;
  int fifthEnd = 0;
  int sixthStart = 0;
  for(i = 0; i < (int)myString.size(); i++){
    if(numSpaces == 4) {
      if(myString[i] == ' ' || myString[i] == '\t') {
        fifthEnd = i;
        sixthStart = i + 1;
        break;
      }
    }
    else if(numSpaces == 3) {
      if(myString[i] == ' ' || myString[i] == '\t') {
        fourthEnd = i;
        fifthStart = i + 1;
        continue;
      }
    }
    else if(numSpaces == 2) {
      if(myString[i] == ' ' || myString[i] == '\t') {
        numSpaces += 1;
        thirdEnd = i;
        fourthStart = i + 1;
        continue;
      }
    }
    else if(numSpaces == 1) {
      if(myString[i] == ' ' || myString[i] == '\t') {
        numSpaces += 1;
        secondEnd = i;
        thirdStart = i + 1;
        continue;
      }
    }
    else if(myString[i] == ' ' || myString[i] == '\t') {
      numSpaces += 1;
      firstEnd = i;
      secondStart = i + 1;
    }
  }
  // convert from string
    int first = atoi(myString.substr(0, firstEnd + 1).c_str());
    int second = atoi(myString.substr(secondStart, secondEnd - secondStart + 1).c_str());
    int third = atoi(myString.substr(thirdStart, thirdEnd - thirdStart + 1).c_str());
    int fourth = atoi(myString.substr(fourthStart, fourthEnd - fourthStart + 1).c_str());
    int fifth = atoi(myString.substr(fifthStart, fifthEnd - fifthStart).c_str());
    int sixth = atoi(myString.substr(sixthStart, myString.size() - sixthStart).c_str());

    if(toReturn == 0){
        return first;
    }
    else if(toReturn == 1){
        return second;
    }
    else if(toReturn == 2){
        return third;
    }
    else if(toReturn == 3){
        return fourth;
    }
    else if(toReturn == 4){
        return fifth;
    }
    else {
      return sixth;
    }
}

int readBenchmark(const char *fileName, routingInst *rst){
    start = high_resolution_clock::now(); //start time used for terminal condition
    ifstream readFile(fileName);
    
    // fill in grid dimentions
    string gridDim;
    getline(readFile, gridDim);
    // find numbers
    int xDim = findXY(gridDim, 1);
    int yDim = findXY(gridDim, 2);

    rst->gx = xDim;
    rst->gy = yDim;

    // fill in capacity
    string capacityLine;
    getline(readFile, capacityLine);
    // capacityLine == capacity XXXX

    // find number
    int cap = findXY(capacityLine, 1);

    rst->cap = cap;

    // fill in num nets
    string numNetsLine;
    getline(readFile, numNetsLine);

    // numNetsLine == num net XXXXXX
    int numNets = findXY(numNetsLine, 2);

    rst->numNets = numNets;

    rst->nets = (net*) malloc(numNets* sizeof(net));
    
    int i;
    // for each net
    for(i = 0; i < numNets; i++) {
    // get num pins
        string numPinsLine;
        getline(readFile, numPinsLine);
        // numPinsLine == netname XXX
        int numPins = findXY(numPinsLine, 1);

        // get all points
        int j;
	      rst->nets[i].pins = (point *) malloc(numPins* sizeof(point));

        for(j = 0; j < numPins; j++) {
          string pointLine;
          getline(readFile, pointLine);
          int x = findXY(pointLine, 0); //pointLine[0] - '0';
          int y = findXY(pointLine, 1); //pointLine[pointLine.length() - 1] - '0';
	        point pin = point();
          pin.x = x;
          pin.y = y;

	        *((*(rst->nets + i)).pins + j) = pin;
        }

        // apply data
	
	
        (*(rst->nets + i)).id = i;
        (*(rst->nets + i)).numPins = numPins;

    }
    



    // calculate edge array
    // naive assume no blockages
    int *numEdges =  (int*)malloc(sizeof(int));
    *(numEdges) = yDim * (xDim - 1) + xDim * (yDim - 1);
    rst->numEdges = *(numEdges);
    rst->edgeCaps = (int*)malloc(*(numEdges) * sizeof(int));
    for(i = 0; i < *(numEdges); i++) {
        *(rst->edgeCaps + i) = cap;
    }

    // find num blockages
    string numBlockagesLine;
    getline(readFile, numBlockagesLine);
    // numBlockagesLine == XXXXX
    int numBlockages = findXY(numBlockagesLine, 0);

    // recalculate edge array for blockages
    for(i = 0; i < numBlockages; i++) {
        string blockageLine;
        getline(readFile, blockageLine);
        int x1 = findXY(blockageLine, 0);
        int y1 = findXY(blockageLine, 1);
        int x2 = findXY(blockageLine, 2);
        int y2 = findXY(blockageLine, 3);
        int new_cap = findXY(blockageLine, 4);
        int edge = getEdgeBetween(x1, y1, x2, y2, xDim, yDim);
        *(rst->edgeCaps + edge) = new_cap;
    }

    readFile.close();


    return 1;
}


int netDecomposition(routingInst *rst) {
  // for each net
  for(int i = 0; i < rst->numNets; i++) {
    
    net *currentNet = (rst->nets + i);
    
    point *oldPinOrdering = (*(rst->nets + i)).pins;
    point *newPinOrdering = (point*) malloc(currentNet->numPins* sizeof(point));
    int numEntered = 0;
    // for each pin in the net
    for(int j = 0; j < currentNet->numPins; j++) {
      if(j == 0) {
        // add first pin as start of the array
        *(newPinOrdering) = *((*(rst->nets + i)).pins);
      }
      else {
        point currentPin = *(newPinOrdering + j - 1);
        // calculate distances between current pin and all others not added
        int distance = 0;
        point closestPin;
        for(int k = 0; k < currentNet->numPins; k++) {
          // skip if old ordering[k] is currentPin or if already in newOrdering
          point checkPin = *(oldPinOrdering + k);
          if(checkPin.x == currentPin.x && checkPin.y == currentPin.y) {
            continue;
          }
          else {
            int notEntered = 1;
            for(int l = 0; l < numEntered; l++) {
              point pinInNewOrdering = *(newPinOrdering + l);
              if(pinInNewOrdering.x == checkPin.x && pinInNewOrdering.y == checkPin.y) {
                notEntered = 0;
                break;
              }
            }
            if(notEntered == 0) {
              continue;
            }
            else {
              // calculate distance
              int currentDistance = sqrt(pow(currentPin.x - checkPin.x, 2) + pow(currentPin.y - checkPin.y, 2));
              if(distance == 0) {
                distance = currentDistance;
                closestPin = checkPin;
              }
              else if(currentDistance < distance) {
                distance = currentDistance;
                closestPin = checkPin;
              }
            }
          }
        }
        *(newPinOrdering + j) = closestPin;
        numEntered += 1;
      }
    }
    free(rst->nets[i].pins);
    rst->nets[i].pins = newPinOrdering;
  }
  return 1;
}

int getEID(point p1, point p2, int gx, int gy) {
    int id=0;
    if (p1.y == p2.y) {
        if (p1.x < p2.x) {
            id = p1.y * (gx-1) + p1.x;
        }
        else {
            id = p2.y * (gx-1) + p2.x;
        }
    }
    else if (p1.x == p2.x) {
        if (p1.y < p2.y) {
            id = p1.y * gx + p1.x + gy * (gx-1);
        }
        else {
            id = p2.y * gx + p2.x + gy * (gx-1);
        }
    }
    return id;
}

point* getEPointer(int id, int gx, int gy){
    point * points = (point *) malloc(2*sizeof(point));
    if(id >= gy * (gx-1)) {
        id -= gy * (gx-1);
        int x = id % gx;
        int y = (id-x)/gx;
        (*(points)).x = x;
        (*(points)).y = y;
        (*(points+1)).x = x;
        (*(points+1)).y = y+1;
    }
    else {
        int x = id % (gx-1);
        int y = (id-x)/(gx-1);
        (*(points)).x = x;
        (*(points)).y = y;
        (*(points+1)).x = x+1;
        (*(points+1)).y = y;
    }

    return points;
}

int solveRouting(routingInst *rst){
  if (rst == NULL) { // if no routingInst passed
    return 0;
  }

  if (rst->gx == 0 || rst->gy == 0 || rst->cap == 0 || rst->numNets == 0 || 
    rst->nets == NULL || rst->numEdges == 0 || rst->edgeCaps == NULL) {
    return 0; //conditions where the routing instant was not properly setup
  }
    
  vector<int*> edgeu;
  for (int i = 0; i < rst-> numNets; i++) { //iterate through nets
    net *opNet = (rst->nets+i);
	  (rst->nets+i)->nroute.numSegs= 0;
	
    (rst->nets+i)->nroute.segments = (segment*) malloc((opNet->numPins-1) * sizeof(segment));//holds segment

    for (int j = 0; j < opNet->numPins-1; j++) { //go through pins of net
      point p1 = *(opNet->pins+j);
      point p2 = *(opNet->pins+j+1); //grab 2 adjacent pins from net
            
      segment newSeg = segment();
      newSeg.p1 = p1;
      newSeg.p2 = p2;
      int *segSize = (int*)malloc(sizeof(int));
            
      vector<int*> combedges;
            
      //if xs are different, go acrosse horizontally
	    point temp1 = p1;
      point temp2 = p1;
      if (p1.x != p2.x) {
        if (p1.x < p2.x) {
          *segSize = p2.x - p1.x;

          temp2.x += 1;
          for(int k = 0; k < *segSize; k++) { //go through edge by edge
            int *edge = (int *) malloc(sizeof(int));
			      *edge = getEID(temp1, temp2, rst->gx, rst->gy); //get id of each edge until p2
            edgeu.push_back(edge);
			      combedges.push_back(edge);
            temp1.x += 1; //move to next edge
            temp2.x += 1;
          }
        }
        else {
          *segSize = p1.x - p2.x;
          temp2.x -= 1;
          for(int k = *segSize; k > 0 ; k--) { //go through edge by edge
              int *edge = (int *) malloc(sizeof(int));
			        *edge = getEID(temp1, temp2, rst->gx, rst->gy); //get id of each edge until p2
              edgeu.push_back(edge);
			        combedges.push_back(edge);
              temp1.x -= 1; //move to next edge
              temp2.x -= 1;
          }
        }
      }
	    temp2.x = temp1.x;
      if (p1.y != p2.y) { //different ys, go upward with edges
		    if (p1.y < p2.y) {
          int iter = p2.y - p1.y;
          temp2.y += 1;
          for(int k = 0; k < iter; k++) { //go through edge by edge
            int *edge = (int *) malloc(sizeof(int));
			      *edge = getEID(temp1, temp2, rst->gx, rst->gy); //get id of each edge until p2
            edgeu.push_back(edge);
			      combedges.push_back(edge);
            temp1.y += 1; //move to next edge
            temp2.y += 1;
          }
        }
        else {
          int iter = p1.y - p2.y;
          *segSize= iter;
          temp2.y -= 1;
          for(int k = iter; k > 0 ; k--) { //go through edge by edge
            int *edge = (int *) malloc(sizeof(int));
			      *edge = getEID(temp1, temp2, rst->gx, rst->gy); //get id of each edge until p2
            edgeu.push_back(edge);
			      combedges.push_back(edge);
            temp1.y -= 1; //move to next edge
            temp2.y -= 1;
          }
        }
      }
      *segSize = combedges.size(); //set number of edges to be size of vector
	    newSeg.numEdges = *segSize;
      newSeg.edges = (int*) malloc(*segSize * sizeof(int));
	    for(int k = 0; k < *segSize; k++) {
        *(newSeg.edges+k) = **(&combedges[0]+k);
	    }
      free(segSize);
            
      *((rst->nets+i)->nroute.segments + j) = newSeg;//add the segment to the route
      (rst->nets+i)->nroute.numSegs= (rst->nets+i)->nroute.numSegs+1;

    }
        
        
  }
  
  //cout<< rst->numEdges << endl;
  rst->edgeUtils = (int *) malloc(rst->numEdges * sizeof(int));
  
  //initialize to all edges used 0 times
  for (int k = 0; k < rst->numEdges; k ++) {
      *(rst->edgeUtils+k) = 0;
  }
  
  //increment edge based on number of edges
  for (int k = 0; k < (int)edgeu.size(); k++) {
    *(rst->edgeUtils + *(edgeu[k])) += 1;
  }
  return 1;
}

void getWeight(int edge, int * weightarray, int* utilhist, int *histarray, routingInst *rst) {
    int util = *(utilhist + edge) - *(rst->edgeCaps + edge); // gets whether usage exceeds cap
    int overflow = max(util, 0); //if util pos, edgecap overflowed
    if (overflow > 0) {
        *(histarray+edge) += 1; //increment history if overflow
    }
    int history = *(histarray+edge);
    int weight = overflow * history; //calculation for edge weight
    *(weightarray + edge) = weight; //update weight in weight array

}

int RRR(routingInst *rst) {
    //initialize weight, history and utilization arrays for edge weights
    int * weightarray = (int *) malloc(rst->numEdges * sizeof(int)); 
    int * histarray = (int *) malloc(rst->numEdges * sizeof(int)); 
    int * utilhist = (int *) malloc(rst->numEdges * sizeof(int)); 
    //initialize hist array to 1 for all edges
    for (int i = 0; i < rst->numEdges; i ++) {
        *(histarray + i) = 1; 
    }
    
    //initialize utilhistory to current edge util
    for (int i = 0; i < rst->numEdges; i ++) {
        *(utilhist + i) = *(rst->edgeUtils+i); 
    }
    
    //initialize all weights in weight array
    for (int i = 0; i < rst->numEdges; i ++) {
        getWeight(i, weightarray, utilhist, histarray, rst); 
    }
    
    //net order calculation 
    vector< pair<int,int>> order; //holds order of nets
    //first int is cost, second is net offset
    for (int i = 0; i < rst->numNets; i++) {
        net * currnet = rst->nets + i;
        route nroute = currnet->nroute;
        int cost = 0;
        for (int j = 0; j < nroute.numSegs; j ++) { //for all segments in the net
            segment * currseg = nroute.segments + j;
            
            for(int k = 0; k < currseg->numEdges; k++) {
                int edge = *(currseg->edges + k);
                cost += *(weightarray + edge);
            }
        }
        order.push_back(make_pair(cost, i)); //adds cost and net pair to order vector
        
    }
    
    //now sort vector descending by cost
    sort(order.rbegin(), order.rend());
    
    //REMOVE: print out of order
    //for(int i = 0; i < (int)order.size(); i++) {
        //cout << order[i].first << ", " << order[i].second << endl;
    //}
    
    //updates the clock before RR
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<minutes>(stop - start);
    /*//Start RRR loop
    //REMOVE Placeholder loop for RRR, will be a for loop instead that breaks if duration.count() > 14
    while(duration.count() < 15) { //while less than 15 minutes
        stop = high_resolution_clock::now();
        duration = duration_cast<minutes>(stop - start);
        cout << duration.count() << endl;
    }*/

	for (int i = 0; i < (int)order.size(); i++) {
		//Rip up
		stop = high_resolution_clock::now();
        duration = duration_cast<minutes>(stop - start);
        if(duration.count() >= 15) {
			break;
		}	

		//cout<< "dur:" <<duration.count()<<endl;


		net *currnet = rst->nets+order[i].second;
		route nroute = currnet->nroute;
		for (int j = 0; j < nroute.numSegs; j ++) { //for all segments in the net
			segment * currseg = nroute.segments + j;
            
            for(int k = 0; k < currseg->numEdges; k++) {
                int edge = *(currseg->edges + k);
				*(rst->edgeUtils+edge) -= 1;
				getWeight(edge, weightarray, utilhist, histarray, rst);
				*(utilhist+edge) += 1;

            }
		}

		//reroute
		for(int j = 0; j < currnet->numPins-1; j++) {
			routeNetAstar(rst, currnet->id, j, j+1);
		}


		//cout<< "thru" << endl;


		for (int j = 0; j < nroute.numSegs; j ++) { //for all segments in the net
			segment * currseg = nroute.segments + j;
            
            for(int k = 0; k < currseg->numEdges; k++) {
                int edge = *(currseg->edges + k);
				*(rst->edgeUtils+edge) += 1;
				getWeight(edge, weightarray, utilhist, histarray, rst);
				//update util history
				//printf("utilhist start\n");
				*(utilhist+edge) += 1;
				//printf("utilhist end\n");	
            }
		}
		
	}
    free(weightarray);
	free(utilhist);
	free(histarray);
    return 1;
}

int writeOutput(const char *outRouteFile, routingInst *rst){
  FILE* writeFile = fopen(outRouteFile, "w");

  int i;
  // for each net
  for(i = 0; i < rst->numNets; i++) {
    fprintf(writeFile, "n%d\n", (*(rst->nets + i)).id);
    
    route finalRoute = (*(rst->nets + i)).nroute;
    int j;
    // for each segment within the route
    for(j = 0; j < finalRoute.numSegs; j++) {
      segment currentSeg = *(finalRoute.segments + j);
      for(int k = 0; k < currentSeg.numEdges; k++) {
        int curedge = *(currentSeg.edges+k);
	      point * pair = getEPointer(curedge, rst->gx, rst->gy);
	      point p1 = *pair;
      	point p2 = *(pair+1);
 	      // write (p1)-(p2)
        fprintf(writeFile, "(%d, %d)-(%d, %d)\n", p1.x, p1.y, p2.x, p2.y);
        free(pair);
      }
      
     
    }
    fprintf(writeFile, "!\n");
  }

  return 1;
}




int release(routingInst *rst){ 
  /*********** TO BE FILLED BY YOU **********/
	int i; 
	int j;
  if(rst != NULL) {
    if (rst->nets != NULL){ // if nets is not NULL
      for(i = 0; i < rst->numNets; i++){ //free each net in the array
        if((*(rst->nets + i)).pins != NULL){
          free((*(rst->nets + i)).pins);
        }
        if((*(rst->nets + i)).nroute.segments != NULL){ // free each segment in the net
          for(j = 0; j < (*(rst->nets + i)).nroute.numSegs; j++){ 
            if((*(rst->nets + i)).nroute.segments[j].edges != NULL){
              free((*(rst->nets + i)).nroute.segments[j].edges);
            }
          }
          free((*(rst->nets + i)).nroute.segments);
        }
      }
      free(rst->nets); // free the array of nets
      free(rst->edgeUtils); // free the edgeUtils array
      free(rst->edgeCaps); // free the edgeCaps array
      free(rst); // free the routingInst
    }
  }
  return(1);
} 

// A* Algortihm

int solveRoutingAstar(routingInst *rst) {

   // Initialize the heap
  PointHash::gy = rst->gy;

  // Initialize the hash table
  // PointHash::hashTable = new PointHash[rst->gx * rst->gy];
  for( int i = 0; i < rst->numNets; i++ )
  {
    
   // printf("net %d\n", i);
   // printf("net %d\n", (*(rst->nets + i)).id);
    rst->nets[i].nroute.numSegs = rst->nets[i].numPins - 1;

    rst->nets[i].nroute.segments = 
      (segment*)malloc( (rst->nets[i].numPins -1)* sizeof(segment)) ;

    for( int j = 0; j < rst->nets[i].numPins - 1; j++)
    {
     // printf("segment %d\n", j);
      int ret = routeNetAstar(rst, i, j, j+1 );
      if( ret == EXIT_FAILURE )
      {

        fprintf( stderr, "failed to route net.\n");
        return EXIT_FAILURE;
      }
    }
  }
  
  return EXIT_SUCCESS;
}




 // A* Algorithm
// TO DO 
int routeNetAstar(routingInst *rst,int netInd ,int SpinInd,int TpinInd) {
  // Initialize the heap
  // 
  unordered_map<point, Node*, PointHash> openSet;
  // unordered_map<point, Node*, PointHash> closedSet;
  priority_queue<costpair, vector<costpair>, greater<costpair> > openSet_pq;
  unordered_map<point, Node*, PointHash> closedSet;

  // Target Node
  Node nT;
  nT.loc = rst->nets[netInd].pins[TpinInd];

  
  // Starting Node
  Node* nS = new Node;
  nS->loc = rst->nets[netInd].pins[SpinInd]; 
  nS->parent = NULL;
  nS->distFromS = 0;

  int dx = abs(nS->loc.x - nT.loc.x);
  int dy = abs(nS->loc.y - nT.loc.y);

  
  nS->distToT = dx + dy;
// printf("nS->distToT %d\n", nS->distToT);
  int approxElems = dx * dy;
  openSet.reserve( approxElems );
  closedSet.reserve( approxElems );
  calcF(*nS);
  // printf("nS->f %d\n", nS->f);
  // 
  openSet.insert(make_pair(nS->loc, nS)); // insert the starting node
  openSet_pq.push( make_pair(nS->F, nS)  ); // remove the starting node
  
  while( !openSet_pq.empty() )
  {
    
    auto top = openSet_pq.top();
	auto current = top.second;
    openSet_pq.pop();

    
    auto it = openSet.find( current->loc );
    if( it == openSet.end() )
    {
      
      continue; 
    }


  //   if( closedSet.find( current->loc ) != closedSet.end() )
  //   {
  //     continue;
  //   }
    if( current->loc.x == nT.loc.x &&  current->loc.y == nT.loc.y)
    {
      int ret = retrace(rst, nS, current, netInd, SpinInd);

  //   printf("ret %d\n", ret);
      deleteMap( closedSet );
      deleteMap( openSet );

      return ret; // success
    }


   // printf("current->loc %d %d\n", current->loc.x, current->loc.y);
    closedSet.insert( std::make_pair(current->loc, current) );
    openSet.erase(it);
    
    
    vector<Node> neighbors;
    
    getNeighbors(rst, neighbors, current, nT );

    
    for( auto neighbor = neighbors.begin(); neighbor != neighbors.end();
              neighbor++ )
    {
      
      auto elem = closedSet.find( neighbor->loc );
      if( elem != closedSet.end() )
      {
        
        continue;
      }
      elem = openSet.find( neighbor->loc );
      if( elem != openSet.end() )
      {
		if(neighbor->F < elem->second->F) {
        
		Node *tmpN = new Node(*neighbor);
		openSet_pq.push( make_pair(tmpN->F, tmpN) );
		openSet.insert( make_pair(tmpN->loc, tmpN));
		}
      }
	  else {
      
		Node *tmpN = new Node(*neighbor);
		openSet_pq.push( make_pair(tmpN->F, tmpN) );
		openSet.insert( make_pair(tmpN->loc, tmpN));
	  } 
    }
  }
  printf( "couldn't find route for netInd: %d\tSpinInd: %d\tTpinInd: %d\n",
      netInd, SpinInd, TpinInd );
  return EXIT_FAILURE; // failure
}





// TO DO
// get the neighbors of the current node
// and put them in the neighbors vector
// return the number of neighbors
// 
int retrace(routingInst *rst, Node* nS, Node* current, int netInd, int segInd) {
Node* tmp = current;

 // 
  rst->nets[netInd].nroute.segments[segInd].p2 = current->loc;
  rst->nets[netInd].nroute.segments[segInd].p1 = nS->loc;
  
  rst->nets[netInd].nroute.segments[segInd].numEdges = current->distFromS;
  rst->nets[netInd].nroute.segments[segInd].edges = 
    (int*)malloc(sizeof(int) * current->distFromS);

  if(rst->nets[netInd].nroute.segments[segInd].edges == NULL)
  {
    fprintf(stderr, "Couldn't malloc for segments[segInd].edges pointer.\n");
    return EXIT_FAILURE;
  }

  int i  = 0;
  double costSum = 0;
  while( tmp->parent != NULL )
  {
    int edgeID = getEID( tmp->loc, (*tmp->parent).loc, rst-> gx, rst->gy );
    rst->nets[netInd].nroute.segments[segInd].edges[i] = edgeID;
    costSum += tmp->G - 1;
  // printf("tmp->loc %d %d\n", tmp->loc.x, tmp->loc.y);
    rst->edgeUtils[edgeID] += 1;

  
    tmp =  tmp->parent;
    i++;
  }

  prevAvgCost = (prevAvgCost + costSum / (double)i ) / 2.0;
  return EXIT_SUCCESS;  
}

// TO DO
void calcF( Node& n ) 
{
  // printf("n.distToT %d\n", n.distToT);
  n.F = n.distFromS + n.distToT;

}

int m_dist( const Node& n1, const Node& n2) { 
  
  return abs(n1.loc.x - n2.loc.x) + abs(n1.loc.y - n2.loc.y);

}

void getNeighbors(routingInst *rst, vector<Node>& neighbors,
                    Node* current, Node T)
{
  if(current->loc.x -1 >= 0) {
    Node newn = Node();
   newn.parent = current; 
    newn.loc.x = current->loc.x -1;
    newn.loc.y = current->loc.y;
    newn.distFromS = current->distFromS + 1;
    newn.distToT = m_dist(newn, T);
    newn.edgeID = getEID(newn.loc, current->loc, rst->gx, rst->gy);
    newn.edgeUtil = *(rst->edgeUtils+newn.edgeID);
    newn.edgeCap = *(rst->edgeCaps+newn.edgeID);
	calcF(newn);

    neighbors.push_back(newn);
  }

  if(current->loc.x +1 < rst->gx) {
    Node newn = Node();
   newn.parent = current; 
    newn.loc.x = current->loc.x +1;
    newn.loc.y = current->loc.y;
    newn.distFromS = current->distFromS + 1;
    newn.distToT = m_dist(newn, T);
    newn.edgeID = getEID(newn.loc, current->loc, rst->gx, rst->gy);
    newn.edgeUtil = *(rst->edgeUtils+newn.edgeID);
    newn.edgeCap = *(rst->edgeCaps+newn.edgeID);
	calcF(newn);

    neighbors.push_back(newn);
  }

  if(current->loc.y -1 >= 0) {
    Node newn = Node();
   newn.parent = current; 
    newn.loc.x = current->loc.x  ;
    newn.loc.y = current->loc.y-1;
    newn.distFromS = current->distFromS + 1;
    newn.distToT = m_dist(newn, T);
    newn.edgeID = getEID(newn.loc, current->loc, rst->gx, rst->gy);
    newn.edgeUtil = *(rst->edgeUtils+newn.edgeID);
    newn.edgeCap = *(rst->edgeCaps+newn.edgeID);
	calcF(newn);

    neighbors.push_back(newn);
  }

  if(current->loc.y +1 < rst->gy) {
    Node newn = Node();
   newn.parent = current; 
    newn.loc.x = current->loc.x ;
    newn.loc.y = current->loc.y+1;
    newn.distFromS = current->distFromS + 1;
    newn.distToT = m_dist(newn, T);
    newn.edgeID = getEID(newn.loc, current->loc, rst->gx, rst->gy);
    newn.edgeUtil = *(rst->edgeUtils+newn.edgeID);
    newn.edgeCap = *(rst->edgeCaps+newn.edgeID);
	calcF(newn);

    neighbors.push_back(newn);
  }

}

void deleteMap( unordered_map<point, Node*, PointHash>& group )
{
  for (unordered_map<point, Node*, PointHash>::iterator it = group.begin(); it != group.end(); ++it) {
    delete it->second;
  }
} 
