// graph.h <Starter Code>
// Caleb Satvedi
/*
Description: This code is the graph class implementation, which
will act as a data structure for our graph data type. It's
container is an unordered map, which holds vertex's and
other unordered maps. The helper functions include a constructor,
adding vertex and edges, and retrieving information about
vertex's neighbors and vertex's in the graph.
*/
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>


using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
unordered_map<VertexT, unordered_map<VertexT, WeightT>> adjList;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //

  graph() {
  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return adjList.size();
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    // counter variable
    int count = 0;
    // use for each loop to go through each unordrerd map
    // in the whole adjList unordrerd map
    for (auto M : adjList) {
      // increment size accordingly
      count += M.second.size();
    }
    return count;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    // return false if said vertex already exists
    if (adjList.count(v) > 0) {
      return false;
    }
    // add a new vertex with an empty unordered map
    unordered_map < VertexT, WeightT> empty;
    adjList.emplace(v, empty);
    return true;
    }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // if either from or to vertex isn't in adjList, return false
    if (adjList.count(from) == 0 || adjList.count(to) == 0) {
      return false;
    } else if (adjList.count(from) != 0 &&
      adjList[from].count(to) != 0) {
      // from's unordered map does have the to vertex
      // updating value
      adjList[from] [to] = weight;
      return true;
    }
    // adding new vertex to inner unordered map
    adjList[from].emplace(to, weight);
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    // if from and to vertex doesn't exist, return false
    if (adjList.count(from) == 0 || adjList.count(to) == 0) {
      return false;
    }
    // if from doesn't go to to, return false
    if (adjList.at(from).count(to) == 0) {
      return false;
    }
    // return weight
    weight = adjList.at(from).at(to);
    return true;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT>  S;
    // vertex doesn't exist in adjaceny list
    if (adjList.count(v) == 0) {
      return S;
    }
    unordered_map <VertexT, WeightT> checker = adjList.at(v);
    // use for each loop to go through each vertex in checker
    // unordered map and add it to set
    for (auto vert : checker) {
      S.insert(vert.first);
    }
    return S;
    }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector <VertexT> v;
    // use foreach loop to add each vertex to vector
    for (auto vert : adjList) {
      v.push_back(vert.first);
    }
    return v;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    // output << "***************************************************" << endl;
    // output << "********************* GRAPH ***********************" << endl;

    // output << "**Num vertices: " << this->NumVertices() << endl;
    // output << "**Num edges: " << this->NumEdges() << endl;

    // output << endl;
    // output << "**Vertices:" << endl;
    // for (int i = 0; i < this->NumVertices(); ++i) {
    //   output << " " << i << ". " << this->Vertices[i] << endl;
    // }

    // output << endl;
    // output << "**Edges:" << endl;
    // for (int row = 0; row < this->NumVertices(); ++row) {
    //   output << " row " << row << ": ";

    //   for (int col = 0; col < this->NumVertices(); ++col) {
    //     if (this->AdjMatrix[row][col].EdgeExists == false) {
    //       output << "F ";
    //     } else {
    //       output << "(T,"
    //         << this->AdjMatrix[row][col].Weight
    //         << ") ";
    //     }
    //   }
    //   output << endl;
    // }
    // output << "**************************************************" << endl;
  }
};
