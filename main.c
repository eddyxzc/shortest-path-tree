#include <stdio.h>
#include<strings.h>
#include<stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
//indexed min heap based priority queue
typedef struct IndexMinPriorityQueue
{
    unsigned int N;// how many elements currently in queue
    unsigned int maxN;// max capacity
    int* pq;
    int* qp;
    int* keys;
}PQueue;

//PQueue constructor
PQueue* PQ_ctor(unsigned int maxCapacity)
{
    if (maxCapacity==0)
    {
        perror("try to creat Pqueue with zero capacity");
        return NULL;
    }
    PQueue* pqueue = malloc(sizeof(PQueue));
    if (pqueue == NULL) {
        perror("PQ_ctor: malloc failed for pqueue\n");
        return NULL;
    }
    pqueue->N = 0;
    pqueue->maxN = maxCapacity;

    pqueue->pq = malloc(sizeof(int) *(maxCapacity+ 1));// array based heap index strat from 1,
    // so aplly one more memry
    if (pqueue->pq == NULL) {
        perror("PQ_ctor: malloc failed for pqueue->pq\n");
        return NULL;
    }

    pqueue->qp = malloc(sizeof(int) * (maxCapacity + 1));
    if (pqueue->qp == NULL) {
        perror("PQ_ctor: malloc failed for pqueue->qp\n");
        return NULL;
    }

    for (unsigned int i = 0; i < maxCapacity + 1; i++)
    {
        pqueue->qp[i]=-1 ;
    }
    pqueue->keys = malloc(sizeof(int) * (maxCapacity + 1));// keys[i] is the object bindding with index i
    if (pqueue->keys == NULL) {
        perror("PQ_ctor: malloc failed for pqueue->keys\n");
        return NULL;
    }
    return pqueue;
}

//PQueue destructor
void PQ_dtor(PQueue* pqueue) {
    if (pqueue != NULL) {
        free(pqueue->pq);
        free(pqueue->qp);
        free(pqueue->keys);
        free(pqueue);
        pqueue = NULL;
    }
}




// comapare the keys mapping to index a and index b
bool greater(PQueue* this, unsigned int a,unsigned int b)
{
    if (this == NULL)
    {
        perror("PQueue is NULL when calling greater() ");
        return false;
    }
    if (a > this->maxN || b > this->maxN) {
        perror("index overflow when calling greater() ");
        return false;
    }
    return this->keys[this->pq[a]] > this->keys[this->pq[b]];
}


// test if pqueue has a key mapped with index k
// i.e. index k has pointted to some key
bool PQ_contains(PQueue* this, unsigned int k)
{
    if (this ==NULL)
    {
        perror("PQueue is NULL when calling contains() ");
        return false;
    }
    return this->qp[k] != -1;
}

// swap two element in heap
void swap(PQueue* this, unsigned int i, unsigned int j) {
    if (this == NULL)
    {
        perror("PQueue is NULL when calling swap() ");
        return ;
    }
    if (i>this->maxN||j>this->maxN) {
        perror("index overflow when calling swap() ");
        return;
    }
     int temp = this->pq[i];
    this->pq[i] = this->pq[j];
    this->pq[j] = temp;

    //update qp too
    this->qp[this->pq[i]] = (int)i;
    this->qp[this->pq[j]] = (int)j;
}

// is parent node is bigger than its child nodes, parent node sink down
/*to move up k/2; to move down k*2(Left Subnode) or k*2+1(Right subnode)*/
void sink(PQueue* this, unsigned int k)
{
    if (this == NULL)
    {
        perror("PQueue is NULL when calling sink() ");
        return;
    }
    if (k > this->maxN) {
        perror("index overflow when calling sink() ");
        return;
    }

    while (2 * k <= this->N) {
        size_t j = 2 * k;

        if (j < this->N && greater(this, j, j + 1))
        {
            // j< this->N means k node has two child
            // // if j=N here, this node k only have left child
            // then we find which child is smaller
            j++;
        }
        // now compare the smaller child with parent
        // if parent is bigger, swap parent with this child node

        if (greater(this, k,j))
        {
            swap(this, k, j);
            k = j;
            // keep sinking
        }
        else {
            // otherwise stop sinking
            break;// parent is smaller than its childs
        }
    }
}

// if child node is bigger than parent, child node swim up
void swim(PQueue* this,unsigned int k)
{
    if (this == NULL)
    {
        perror("PQueue is NULL when calling swim() ");
        return;
    }
    if (k > this->maxN) {
        perror("index overflow when calling swim() ");
        return;
    }
    // k/2 is the parent node
    // if k==1, this node is already the root node
    while (k > 1 && greater(this, k / 2, k))
    {
        //parent node is bigger, swim to parent node
        swap(this, k / 2, k);
        k = k / 2;
    }
}

// insert a key bindding it to index k
void PQ_insert(PQueue* this, unsigned int k, int key)
{
    if (this == NULL)
    {
        perror("PQueue is NULL when calling PQ_insert() ");
        return;
    }
    if (k > this->maxN) {
        perror("index overflow when calling PQ_insert() ");
        return;
    }
    if (PQ_contains(this, k)) {
        //already have this index k bindding to some key
        return;
    }
    else {
        // put the new key at the bottom of the heap
        this->N++;
        this->pq[this->N] = (int)k;
        this->qp[k] = (int)this->N;
        this->keys[k] = key;
        // then swim up this new key
        swim(this, this->N);
    }
}

// return the index at the front
int PQ_peekIndex(PQueue* this)
{

    return this->pq[1];
}

// return the key at the front
int PQ_peekKey(PQueue* this)
{

    return this->keys[PQ_peekIndex(this)];
}

// pop the front, the minimum key
int PQ_pop(PQueue* this)
{

    if (this == NULL)
    {
        perror("PQueue is NULL when calling PQ_pop()\n ");
        return -1;
    }
    if (this->N == 0) {
        perror("pop from empty queue");
        return -1;
    }

    int popedkey = PQ_peekKey(this);
    unsigned int popedindex = PQ_peekIndex(this);


    swap(this, 1, this->N);// let the last element. to the heap head.
    // then sink this element, rehipy the heap
    this->N--;
    this->qp[popedindex] = -1;// unmapping this index, as it not valid anymore
    sink(this,1);
    return popedkey;

}

// replace the key bidding to index k
void PQ_replace(PQueue* this, unsigned int k, int key)
{
    if (this == NULL)
    {
        perror("PQueue is NULL when calling PQ_replace()\n ");
        return;
    }

    if (!PQ_contains(this, k)) {
        perror("PQ_replace: no such k\n");
        return;
    }
    this->keys[k] = key;
    swim(this, this->qp[k]);
    sink(this, this->qp[k]);
    return;
}

typedef  int WeightT;
//struct WeightedDirectedEdgeType
// from v--->w
typedef struct WeightedDirectedEdgeType
{
    size_t v;// vertex v, from v---->w,. source vertex
    size_t w;// vertex w, endpoint
    WeightT weight;
}WDEdge;

/*adjacencylist of A vertex, contains all edges from this vertex */
// like C++ vector<WDEdges>
typedef struct WeightedDirectedEdgeAdjacencyListType
{
    WDEdge* edges;
    size_t space;
    size_t num;

}WDAjList;

typedef struct WeightedDirectedGraph
{
    size_t v;// number of vertex
    size_t e;// number of egeds
    WDAjList** adjlist; // adjacencylist number adjlist[v] is the vertex v's adjacency list
}WDgraph;


WDAjList* WDAjList_ctor()
{
    WDAjList* this = malloc(sizeof(WDAjList));
    if (this == NULL)
    {
        perror("WDAjList_ctor: allocate memory for this failed!\n");
        return NULL;
    }
    this->edges = malloc(sizeof(WDEdge)*4); // preallocate 4 space
    if (this->edges == NULL)
    {
        perror("WDAjList_ctor: allocate memory for this->edges failed!\n");
        free(this);
        return NULL;
    }
    this->space = 4;
    this->num = 0;
    return this;
}

void WDAjList_dtor(WDAjList* this)
{
    if (this != NULL) {
        free(this->edges);
        free(this);
        this = NULL;
    }
}

// push back an edge, edge will be copyied into WDAjlist
// assign index of this new element
// return 0 succeed or -1 for error
int WDAjList_pushback(WDAjList* this, WDEdge* edge, size_t* index) {
    if (this == NULL) {
        printf("NULL pointer in WDAjList_pushback");
        return -1;
    }

    if (this->num < this->space)
    {
        //no over flow
    }
    else {
        // overflow, extend space first, two times large now
        WDEdge* newedgespace = realloc(this->edges,sizeof(WDEdge) * this->space * 2);
        if (newedgespace == NULL)
        {
            perror(" extend memory failed in WDAjList_pushback!\n");
            return -1;
        }

        // mount new space
        this->edges = newedgespace;
        this->space *= 2;
    }


    *index = this->num;
    this->edges[*index].v = edge->v;
    this->edges[*index].w = edge->w;
    this->edges[*index].weight = edge->weight;
    this->num++; // next insert element will be here
    return 0;
}

// constructor of WDgraph
WDgraph* WDgraph_ctor(size_t numVertex)
{
    if (numVertex == 0) {
        perror("WDgraph_ctor: Create WDgraph with Zero vertex\n ");
        return NULL;
    }
    WDgraph* this = malloc(sizeof(WDgraph));
    if (this == NULL) {
        perror("Allocate graph failed! in WDgraph_ctor\n");
        return NULL;
    }
    this->v = numVertex;
    this->e = 0;
    this->adjlist = malloc(sizeof(WDAjList*) * numVertex);
    if (this->adjlist == NULL) {
        perror("Allocate adjlist[] failed! in WDgraph_ctor\n");
        free(this);
        return NULL;
    }
    //allocate memory for each vertex's adjancncy list
    for (size_t i = 0; i < this->v; i++)
    {
        this->adjlist[i]=WDAjList_ctor();
        if (this->adjlist[i] == NULL) {
            perror("Allocate adjlist[i] failed! in WDgraph_ctor\n");
            free(this);
            return NULL;
        }
    }
    return this;
}

//destructor of WDgraph
void WDgraph_dtor(WDgraph* this) {
    if (this == NULL)return;
    for(size_t i=0;i<this->v;i++)
    { WDAjList_dtor(this->adjlist[i]);}
    free(this->adjlist);
    free(this);
}

// simple add weighted directed edge to WDgraph v---->w, weight
// no check whether this edge is in the graph already.
// failed return -1
int WDgraph_addEdge(WDgraph* this, unsigned int v, unsigned int w, WeightT weight)
{
    if (this == NULL)
    {
        perror("Calling WDgraph_addEdge with NULL pointer\n");
        return -1;
    }
    if (weight < 0) {
        printf("Negative weight detected!\n");
    }
    WDEdge newedge;
    newedge.v = v;
    newedge.w = w;
    newedge.weight = weight;
    size_t index;
    if (WDAjList_pushback(this->adjlist[v],&newedge,&index) != 0) {
        perror("WDAjList_pushback failed in WDgraph_addEdge\n");
        return -1;
    }
    this->e++;// |edge| ++
    return 0;
}

typedef struct  ShortestPathTree_AUX_DATA_TYPE
{
    WDEdge* edges;//  edges[v] is a edge in Shortest Path Tree
    WeightT* dists;// dists[v] is the distance from source to vertex v,
    PQueue* pqueue;// min indexed proirity queue

    unsigned int sourceVertex;// source vertex
    WDgraph* graph;// refer to an OUTSIDE directed Weighted graph, NO new or delete
}SPTdata;

//Construct SPTdata from a vbalid non negtive weighted directed graph, and a source vertex
//don't free graph before free SPTdata
SPTdata* SPTdata_ctor(WDgraph* graph, unsigned int sourceVertex) {
    if (graph == NULL) {
        perror(" Creat SPTdata from a NULL graph pointer\n");
        return NULL;
    }

    SPTdata* this = malloc(sizeof(SPTdata));
    if (this == NULL) {
        perror("Allocate memory failed in SPTdata_ctor\n ");
        return NULL;
    }
    this->edges = malloc(sizeof(WDEdge) * graph->v);// at most All vertex on Shortest Path Tree
    if (this->edges == NULL) {
        perror("Allocate memory failed in SPTdata_ctor,this->edges == NULL\n ");
        return NULL;
    }

    this->dists = malloc(sizeof(WeightT) * graph->v);// all vertex to source's distances
    if (this->dists == NULL) {
        perror("Allocate memory failed in SPTdata_ctor this->dists==NULL\n ");
        return NULL;
    }
    for (size_t i = 0; i < graph->v; i++)
    {
        this->dists[i] = INT_MAX; // set all vertex to source as infinity
        this->edges[i].weight= UINT_MAX;//
    }

    this->dists[0] = 0;//source too source is 0
    this->pqueue = PQ_ctor(graph->v);
    this->graph = graph;
    this->sourceVertex=sourceVertex;
    return this;
}


//SPTdata destructor, note: graph won't be free
void SPTdata_dtor(SPTdata* this)
{
    if (this != NULL)
    {
        free(this->edges);
        free(this->dists);
        PQ_dtor(this->pqueue);
        free(this);
        this = NULL;
    }
}

void SPT_relaxEdge(SPTdata* this, int vertex)
{
    const WDAjList* adjlist = this->graph->adjlist[vertex];
    PQueue* pqueue = this->pqueue;
    for (size_t i = 0; i < adjlist->num; i++)
    {
        unsigned int w = adjlist->edges[i].w;
        if (this->dists[w] > this->dists[vertex] + adjlist->edges[i].weight)
        {
            this->dists[w] = this->dists[vertex] + adjlist->edges[i].weight;
            this->edges[w] = adjlist->edges[i];
            if (PQ_contains(pqueue, w))
            {
                PQ_replace(pqueue, w, this->dists[w]);
            }
            else {
                PQ_insert(pqueue, w, this->dists[w]);

            }
        }

    }

}


int CalculateSPT(SPTdata* sptdata)
{
    PQueue* pqueue = sptdata->pqueue;
    // first put source into the queue
    PQ_insert(pqueue, sptdata->sourceVertex, 0);
    //then relax all edges connected to source, directly and un directly
    while (pqueue->N != 0)
    {
        int nextVertex = PQ_peekIndex(pqueue);
        PQ_pop(pqueue);
        SPT_relaxEdge(sptdata, nextVertex);//releax edges from nextVertex
    }

    return 0;
}
/*
     find the longest path inside SPTree
	 the weight is the answer to this problem
	*/
int get_longest_path(SPTdata* this)
{
    int max=-1;
    for (int i = 0; i < this->graph->v; i++) {
        max=max>this->dists[i]?max:this->dists[i];
    }
    return max;
}


// max number of chars at one time when parse input from stdin
#define LINEMAX 1024

//get a line from file, return -1 if failed
// lineptr is the line string, it is allocated inside but need to be free outside
// n is the number of chars in string(not includ the last "\0")
int mygetline(char** lineptr, size_t* n, FILE* stream)
{
    char* buffer = *lineptr;
    char* relocate_buffer = buffer;

    buffer = malloc(LINEMAX+1);

    if (buffer == NULL)
    {
        printf("allocate memory failed in getline!\n");
        return -1;
    }
    size_t linemax = LINEMAX + 1;
    size_t numread = 0;
    char c;
    while (numread< linemax)
    {
        c = (char)fgetc(stdin);
        if (c == EOF)break;// end of file
        buffer[numread] = c;
        numread++;
        if (c == '\n')// end of line
            break;
        if (numread == linemax-1) {
            //try to extend memory to 2 times
            relocate_buffer = realloc(buffer, 2 * linemax);
            if (relocate_buffer == NULL) {
                free(buffer);
                printf("Error when extend memory in getline!\n");
                return -1;
            }else {
                buffer = relocate_buffer;
                linemax *= 2;
            }
        }

    }
    buffer[numread] = '\0';
    *lineptr = buffer;
    *n = numread;
    return 0;
}

//parse adjacency matrix,read a line string, get the a row of matrix, return the number of column
// return -1 if failed. if the char is 'x' the value in matrix will be -1
// the matrx_row_ptr must be already allocated;
int ParseWeight(char* lineptr, int* matrx_row) {
    size_t offset = 0;// line str parse offset
    size_t conlumn = 0; //conlumn index of this row
    while (lineptr[offset] != '\n')
    {
        if (lineptr[offset] == ' ')
        {
            offset++;
            continue;
        }

        if (lineptr[offset] == 'x')
        {
            offset++; // update str offset
            //do something
            // /*Debug*/printf("%c\n", 'x');
            matrx_row[conlumn] = -1;
            conlumn++;// update conlumn index
            continue;
        }
        char* nextchar;
        int weight = (int)strtol(&lineptr[offset], &nextchar, 10);// try to get the number
        if (weight == 0) {
            perror("ParseWeight failed!\n");
            return -1; }// farse failed
        offset = nextchar - lineptr;// update offset

        ///*Debug*/printf("%d\n", weight);
        matrx_row[conlumn] = weight;
        conlumn++;// update conlumn index
    }
    return conlumn;
}
int main(int argc, char** argv)
{


    char* lineptr = NULL;
    size_t linelength = 0;
    //read number of cities
    int rt = -1;
    rt=mygetline(&lineptr, &linelength, stdin);
    if (rt == -1)
    {
        printf("Read number of city failed!");
        return -1;
    }
    char* nextchar;
    int ncity = (int)strtol(lineptr, &nextchar, 10);// try to get the number
    free(lineptr);
    if( ncity < 1) {
        printf("INPUT ERROR, please input city number N first\n");
        return -1;
    }

    // now read the adjacency matrix(should be n-1 lines)

    //adjacency matrix
    int** adj_matrix = malloc(sizeof(int*) * (ncity - 1));
    if (adj_matrix == NULL) {
        perror("Adjacency matrix allocation failed!\n");
        return -1;
    }

    for (int i = 0; i < ncity-1; i++)
    {
        //printf("int i: %d ",i);
        rt=mygetline(&lineptr, &linelength, stdin);
        if (rt == -1)
        {
            printf("Read adjacency lines failed!");
            return -1;
        }
        // allocate memory for adj_matrix[i]
        adj_matrix[i] = malloc(sizeof(int) * (i+1));
        if (adj_matrix[i] == NULL) {
            printf("Adjacency matrix row %d allocation failed!\n", i);
            return -1;
        }

        rt=ParseWeight(lineptr, adj_matrix[i]);
        if (rt==-1)
        {
            printf("Adjacency matrix row %d ParseWeight failed!\n", i);
            return -1;
        }
        free(lineptr);//lineptr allocated inside getline
    }

    // construct a graph

    WDgraph* graph = WDgraph_ctor(ncity);
    if (graph == NULL) {
        perror("Allocate graph failed!\n");
        return -1;
    }
    // add edges
    for (int irow = 0; irow < ncity-1; irow++)
    {
        for (int jcol = 0; jcol <= irow; jcol++)
        {
            if(adj_matrix[irow][jcol]==-1)continue;
            if (WDgraph_addEdge(graph, irow + 1, jcol, adj_matrix[irow][jcol]) != 0) {
                perror("Add Edge failed!\n");
                return -1;
            }
            // add w->v too, as we can reach cities from both directions
            if (WDgraph_addEdge(graph, jcol, irow+1, adj_matrix[irow][jcol]) != 0) {
                perror("Add Edge failed!\n");
                return -1;
            }

            ///*debug*/printf("%d--%d-->%d......%d--%d-->%d\n", irow+1, adj_matrix[irow][jcol], jcol, jcol,
             //               adj_matrix[irow][jcol], irow+1);
        }
    }

    // construct shortest path tree
    SPTdata* sptree = SPTdata_ctor(graph, 0);
    if (sptree == NULL) {
        perror("Allocate SPTdata failed!\n");
        return -1;
    }
    // canculate shortest path tree
    CalculateSPT(sptree);


    //get the answer
    printf("%d\n",get_longest_path(sptree));

    // free resources and return
    SPTdata_dtor(sptree);
    WDgraph_dtor(graph);
    for (int i = 0; i < ncity-1; i++)
    {
        /*DEBUG
        for (int j = 0; j <= i; j++)
        {
            printf("%d ", adj_matrix[i][j]);
        }*/
        free(adj_matrix[i]);
        //printf("\n");
    }
    free(adj_matrix);
    return 0;
}
