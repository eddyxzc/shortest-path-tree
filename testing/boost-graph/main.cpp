#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

using namespace boost;


#define LINEMAX 1024
//get a line from file, return -1 if failed
// lineptr is the line string, it is allocated inside but need to be free outside
// n is the number of chars in string(not includ the last "\0")
int mygetline(char** lineptr, size_t* n, FILE* stream)
{
    char* buffer = *lineptr;
    char* relocate_buffer = buffer;

    buffer = static_cast<char *>(malloc(LINEMAX + 1));

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
            relocate_buffer = static_cast<char *>(realloc(buffer, 2 * linemax));
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



int main(int, char*[])
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
    int** adj_matrix = static_cast<int **>(malloc(sizeof(int *) * (ncity - 1)));
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
        adj_matrix[i] = static_cast<int *>(malloc(sizeof(int) * (i + 1)));
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








    typedef adjacency_list< listS, vecS, undirectedS, no_property,
            property< edge_weight_t, int > >
            graph_t;
    typedef graph_traits< graph_t >::vertex_descriptor vertex_descriptor;
    typedef std::pair< int, int > Edge;


    std::vector<Edge> edge_array;
    std::vector<int> weights;
    // add edges
    for (int irow = 0; irow < ncity-1; irow++)
    {
        for (int jcol = 0; jcol <= irow; jcol++)
        {
            if(adj_matrix[irow][jcol]==-1)continue;
            edge_array.emplace_back(irow + 1, jcol);
            weights.push_back(adj_matrix[irow][jcol]);

            /*// add w->v too, as we can reach cities from both directions
            edge_array.emplace_back(jcol, irow+1);
            weights.push_back(adj_matrix[irow][jcol]);*/

            ///*debug*/printf("%d--%d-->%d......%d--%d-->%d\n", irow+1, adj_matrix[irow][jcol], jcol, jcol,
            //               adj_matrix[irow][jcol], irow+1);
        }
    }


    int num_arcs = edge_array.size();
    graph_t g(edge_array.data(), edge_array.data()+num_arcs, weights.data(), ncity);
    property_map< graph_t, edge_weight_t >::type weightmap
            = get(edge_weight, g);
    std::vector< vertex_descriptor > p(num_vertices(g));
    std::vector< int > d(num_vertices(g));
    vertex_descriptor s = vertex(0, g);

    /*dijkstra_shortest_paths(g, s,
                            predecessor_map(boost::make_iterator_property_map(
                                    p.begin(), get(boost::vertex_index, g)))
                                    .distance_map(boost::make_iterator_property_map(
                                            d.begin(), get(boost::vertex_index, g))));*/

    dijkstra_shortest_paths(g, s,
                            distance_map(boost::make_iterator_property_map(
                                            d.begin(), get(boost::vertex_index, g))));

    //std::cout << "distances and parents:" << std::endl;
    graph_traits< graph_t >::vertex_iterator vi, vend;
    std::vector<int> dist[100];
    /*for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi)
    {
        std::cout << "distance(" << *vi << ") = " << d[*vi] << ", ";
        std::cout << "parent(" << *vi << ") = " << p[*vi]
                  << std::endl;
    }*/
    auto max=*std::max_element(d.begin(),d.end());
    std::cout << max << std::endl;


    return EXIT_SUCCESS;
}