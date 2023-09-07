// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.


#ifndef ECE556_H
#define ECE556_H

#include <stdio.h>
#include <type_traits>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include <math.h>


 /**
  * A structure to represent a 2D Point. 
  */
 struct point
 {
   int x ; /* x coordinate ( >=0 in the routing grid)*/
   int y ; /* y coordinate ( >=0 in the routing grid)*/
   
	bool operator==(const point& p) const
  {
    return (x == p.x && y == p.y);
  }

 };


  /**
  * A structure to represent a segment
  */
 typedef struct
 {
   point p1 ; 	/* start point of a segment */
   point p2 ; 	/* end point of a segment */
   
   int numEdges ; 	/* number of edges in the segment*/
   int *edges ;  	/* array of edges representing the segment*/
   
 } segment ;
 
 
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
   This function creates a routing solution
   input: pointer to the routing instance
   output: 1 if successful, 0 otherwise (e.g. the data structures are not populated) 
*/
int solveRouting(routingInst *rst);
  
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

int getEID(point p1, point p2, int gx, int gy);

point* getEPoint(int id, int gx, int gy);

    /* void getWeight(int edge, int * weightarray, int *histarray, routingInst *rst);
        Calculates weight of given edge for history based Rip up and Reroute. 
        Edge is the edge that the weight is being calculated for
        weightarray is an array of all the weights of every edge, updated in function
        histarray is an array of all the previous histories of each edge in prior iterations, updates in function
        rst is a routing instance holding edge capacities and utilizations
    */
void getWeight(int edge, int * weightarray, int *histarray, routingInst *rst);

    /* int RRR(routingInst *rst);
        executes route rip and replace
        first gets all weights of every edge, then orders the nets by cost of 
        weights of edges in route
    */
int RRR(routingInst *rst);

    /*int netDecomposition(routingInst *rst);
      reorders pins within each net
      ensures that overall pathing is shorter through initial solution
    */
int netDecomposition(routingInst *rst);



/* parts of namespace used */
using std::vector;
using std::pair;
using std::unordered_map;
using std::priority_queue;


const double MAX_HEUR_COST = 20.0;
const double K_OVF_MULT = 1.0;
const int K_EDGE_CAPS = 1;
const int K_EDGE_UTIL = 1;

const double unitDist = 0.99;


static double prevAvgCost = 0;

// used to store the edges of the grid

namespace{

struct Node
{
  point loc;   // location of node
  Node* parent; // parent node
  int edgeID;    // edgeID of edge used to reach this node
  int distFromS; // distance from source
  int distToT;   // distance to target
  int edgeUtil;  // edge utilization
  int edgeCap;   // edge capacity

  double F;    //total cost of this node (G + H)
  double G;    // cost of getting to this node from source



  Node() : loc(point()), parent(NULL), edgeID(0), distFromS(-1),
            distToT(-1), edgeUtil(0), edgeCap(0), F(-1), G(-1)
  {

	loc.x = -1;
	loc.y = -1;
  }

  Node( const Node& other ): loc(other.loc), parent(other.parent), edgeID(other.edgeID),
                      distFromS(other.distFromS), distToT(other.distToT),
                      edgeUtil(other.edgeUtil), edgeCap(other.edgeCap),
                      F(other.F), G(other.G) {}

  bool operator==(const Node& n) const
  {
    return (loc.x == n.loc.x && loc.y == n.loc.y);
  }
  bool operator!=(const Node&n) const
  {
    return !(loc.x == n.loc.x && loc.y == n.loc.y);
  }
};



class CompareNodeCost
{
  public:
    bool operator()( const Node* n1, const Node* n2 ) const
    {
     
      return n1->F > n2->F;
    }
};      


class PointHash
{
  public:
    std::size_t operator()(const point& pt) const
    {
      return pt.y + gy * pt.x;
    }

    static int gy;
};

int PointHash::gy = 0;

} //end namespace


void calcF( Node& n ) ;
void calcG( Node& n ) ;

void getNeighbors(routingInst *rst, vector<Node>& neighbors,
                    Node* current, Node T);

int retrace(routingInst *rst, Node* nS, Node* current, int netInd, int segInd) ;

int solveRoutingAstar(routingInst *rst);

int routeNetAstar(routingInst *rst,int netInd,int SpinInd,int TpinInd);
int m_dist( const Node& n1, const Node& n2);

void deleteMap( unordered_map<point, Node*, PointHash>& group );

    

#endif // ECE556_H