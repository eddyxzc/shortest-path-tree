# shortest-path-tree



caculate shortest path tree with Dijkstra’s algorithm

to build:
make

to run:
./sssp


to test:(need boost graph library and g++ cmake)
cd testing
./test.sh


The complexity should be O(V+E) as I use indexed min-heap as a priority queue.
OOP technology is used to make dynamic memory more manageable.
It was tested with Valgrind. No memory leak was found.

reference: Algorithms (4th Edition) by Robert Sedgewick , Kevin Wayne 
API：
## data structure

 	1. IndexMinPriorityQueue
      	1. `PQ_ctor` : //Constructor
      	2. `PQ_dtor` : //Destructor
      	3. `Sink` : // is parent node is bigger than its child nodes, parent node sink down
      	4. `Swim` : // if child node is bigger than parent, child node swim up
      	5. `PQ_insert` : // insert a key bindding it to index k
      	6. `PQ_peekIndex` : // return the index at the front
      	7. `PQ_peekKey` : // return the key at the front
      	8. `PQ_pop` : // pop the front, the minimum key
      	9. `PQ_replace` : // replace the key bidding to index k
      	10. `greater` : // comapare the keys mapping to index a and index b
      	11. `PQ_contains `  : // test if pqueue has a key mapped with index k
      	12. `swap  `  : // swap two element in heap
 	2. WeightedDirectedEdgeType
 	3. WeightedDirectedEdgeAdjacencyListType
      	1. `WDAjList_ctor` : //Constructor
      	2. `WDAjList_dtor` : //Destructor
      	3. `WDAjList_pushback` : // push back an edge, edge will be copyied into WDAjlist
 	4. WeightedDirectedGraph
      	1. `WDgraph_ctor` : //Constructor
      	2. `WDgraph_dtor`: //Destructor
      	3. `WDgraph_addEdge` : // simple add weighted directed edge to WDgraph v---->w, weight
 	5. ShortestPathTree_AUX_DATA_TYPE
      	1. `SPTdata_ctor` : //Constructor
      	2. `SPTdata_dtor` : //Destructor
      	3. `SPT_relaxEdge`
      	4. `CalculateSPT`
      	5. `get_longest_path` : // find the longest path inside SPTree, the weight is the answer to this problem



## input

1. `mygetline`

   ```
   //get a line from file, return -1 if failed
   // lineptr is the line string, it is allocated inside but need to be free outside
   // n is the number of chars in string(not includ the last "\0")
   ```

2. `ParseWeight`

   ```
   //parse adjacency matrix,read a line string, get the a row of matrix, return the number of column
   // return -1 if failed. if the char is 'x' the value in matrix will be -1
   // the matrx_row_ptr must be already allocated;
