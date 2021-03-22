#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/graph_traits.hpp>

typedef boost::adjacency_list<> Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;

using namespace boost;
#include <random>
#include <vector>
#include<string>
#include<iostream>
#include<fstream>
using std::vector;
using std::string;
int main()
{
    int nodes=100;
    std::mt19937 weightgenerator { std::random_device{}() };
    std::uniform_int_distribution<int> dist(1, nodes);
    boost::minstd_rand prng(std::random_device{}());
    // Create graph with 100 nodes and edges with probability 0.05
    Graph g(ERGen(prng, nodes, 0.5), ERGen(), nodes);
    auto vertex_idMap = get(boost::vertex_index, g);
    boost::graph_traits <Graph>::vertex_iterator i, end;
    boost::graph_traits <Graph>::adjacency_iterator ai, a_end;


    vector<vector<int>> matrix(nodes);
    for (auto& rows:matrix)
    {
        rows=vector<int>(nodes,-1);
    }


    for (boost::tie(i, end) = vertices(g); i != end; ++i) {
        //std::cout << vertex_idMap[*i] << ": ";

        for (boost::tie(ai, a_end) = adjacent_vertices(*i, g); ai != a_end; ++ai) {
            weightgenerator.seed(std::random_device{}());
            matrix[vertex_idMap[*i]][vertex_idMap[*ai]]=dist(weightgenerator);
            //std::cout << vertex_idMap[*ai]<<'('<<matrix[vertex_idMap[*i]][vertex_idMap[*ai]]<<')';
            //if (boost::next(ai) != a_end)
               // std::cout << ", ";
        }
        //std::cout << std::endl;
    }

    // output matrix

    std::ofstream of("tempMatrix.data",std::ios::trunc);
    of<<std::to_string(nodes)<<std::endl;
    for (int i=1;i<nodes;i++) {
        for(int j=0;j<i;j++){
            if(matrix[i][j]!=-1)
            {
                of<<std::to_string(matrix[i][j]);
            }else{
                of<<'x';
            }
            if((j+1)==i){
                of<<std::endl;
            }else {
                of<<' ';
            }
    }
    }

    return 0;
}