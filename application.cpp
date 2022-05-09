// application.cpp <Starter Code>
/*
Program 7
Caleb Satvedi
This is the code implementation for the open street maps class. In this file, we 
have the code for how to searching for a building, finding the building 
in the middle of 2 other buildings, setting the nodes for said buildings,
formulating a dijkstra algroithm to find the shortest path between 2 nodes,
and finding the path between 2 nodes, as well as functions that call other
functions and those that print statements. Creative component right above
application()
*/
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include <queue>
#include <stack>
#include <algorithm>  // for find()
#include <cmath>  // for absolute value


#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;
const double INF = numeric_limits<double>::max();
const int stepsInMile = 2000;

// We create this class to help with
// the priotity queue in dijstrkas algorithm funciton

class prioritize  // you could also use a struct
{
public:
// the inputs are 2 pairs of type long long and double
  bool operator()(const pair<long long, double>& p1,
  const pair<long long, double>& p2) const  {
    return p1.second > p2.second;
  }
};

// We search for a building that has the same name as the query sent in
BuildingInfo searchBuilding(string query,
vector<BuildingInfo> Buildings) {
BuildingInfo ourBuild;
// use for-each loop to go through buildings vector
for (auto builder : Buildings) {
  // check if abbreviation is same
  if (builder.Abbrev == query) {
    // if so, set that to the building to be returned
    ourBuild = builder;
    break;
  }
  // check if full name is the same
  if (builder.Fullname.find(query) != string::npos) {
    // if so, set that to be returned
    ourBuild = builder;
    break;
  }
}
// return the building
  return ourBuild;
}


// We find the building in the middle of two other buildings
BuildingInfo nearestBuildingSet(BuildingInfo p1, BuildingInfo p2,
vector<BuildingInfo> Buildings, set <string> excluded) {
BuildingInfo closest;
// Find the coordinates which are betweeen the middle of the two buildings
Coordinates middle = centerBetween2Points(p1.Coords.Lat, p1.Coords.Lon, p2.Coords.Lat, p2.Coords.Lon);
  // set min to infinity
  double min = INF;
  // use for each loop to go through each building in buildings vector
  for (auto build : Buildings) {
    // if the building name is in the excluded building, skip this iteration
    if (excluded.count(build.Fullname) != 0) {
      continue;
    }
    // find the distance between the building and middle cooridinate
    double distance = distBetween2Points(middle.Lat, middle.Lon,
      build.Coords.Lat, build.Coords.Lon);
      // if the distance is the smallest, update holders
    if (distance < min) {
      min = distance;
      closest = build;
    }
  }
  // return the closest building
return closest;
}

// This function find the node which is closest to the building we send in
long long nearestNode(BuildingInfo b, vector<FootwayInfo> Footways,
map<long long, Coordinates> Nodes) {
  // set holder variables to closest node to building
  long long closest;
  double min = INF;
  // for each loop through footways vector
  for (auto fw : Footways) {
    // for each loop through nodes vector in single footway struct
    for (auto n : fw.Nodes) {
      // find the distance between the building and the node we are on
      double distance = distBetween2Points(b.Coords.Lat, b.Coords.Lon, Nodes[n].Lat, Nodes[n].Lon);
      // if the distance is smallest so far, we update it
      if (distance < min) {
        min = distance;
        closest = Nodes[n].ID;
      }
    }
  }
  return closest;
}

// This function implements dijkstra's algorighm
// to find the shortest path the start node
// and all the other nodes
void dijkstra(graph<long long, double> G,
map<long long, double>& distances,
map<long long, long long>& predecessors, long long start) {
  // Create a priority queue to list the nodes
  // and their distances from the building
  priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> unvisitedQueue;
  // Use getverticies function in graph class to get all verticies in the graph
  vector <long long> allvertex = G.getVertices();
  // for each loop through each vertex in allvertex vector
  for (auto vert : allvertex) {
    // set each each vertex to respective value for
    // distances, predecessors, and unvisitedQueue
    distances.emplace(vert, INF);
    predecessors.emplace(vert, 0);
    unvisitedQueue.push(make_pair(vert, INF));
      }
    // set start distance to 0
  distances[start] = 0;
  // created a visited set
  set<long long> visited;
  // make sure the unvisitedQueue has 0 distance for start
  unvisitedQueue.push(make_pair(start, 0));
  // while loop for when unvisitedQueue is still full
  while (!unvisitedQueue.empty()) {
    // take off top element from unvisited queue
    long long currV = unvisitedQueue.top().first;
    double topW = unvisitedQueue.top().second;
    unvisitedQueue.pop();
    // if statements for to know when to
    // break, continue, or go on with while loop
    if (topW == INF) {
      break;
    } else if (visited.count(currV) != 0) {
      continue;
    } else {
      // use neighbors function from graph class
      // to get currV's neighbors
      set<long long> adjacent = G.neighbors(currV);
      // go through each neighbor in a for-each loop
      for (auto adj : adjacent) {
        // get edge weigh from currV to adjacent
        double edgeWeight = 0.0;
        G.getWeight(currV, adj, edgeWeight);
        // make edge weight from startV to adj
        double altDistance = distances[currV] + edgeWeight;
        // if the new edgeWeight is less than the current one
        if (altDistance < distances[adj]) {
          // then we update the values in the maps and queue
          distances[adj] = altDistance;
          predecessors[adj] = currV;
          unvisitedQueue.push(make_pair(adj, altDistance));
                  }
                  // add currV to the visited set
        visited.insert(currV);
        }
      }
    }
  }

// this function returns the path
// for the vertex we are trying to go to
vector <long long> getPath(map<long long, long long> predecessors,
long long endVertex) {
  // create a stack for the path
  stack <long long> pathBack;
  long long currV = endVertex;
  // add last vertex to stack and go to next one
  // according to predecessors map
  while (currV != 0) {
    pathBack.push(currV);
    currV = predecessors[currV];
      }
    // create a vector to formulate the path correctly
  vector <long long> path;
  while (!pathBack.empty()) {
    // put the first element of the stack in the back of the vector
    currV = pathBack.top();
    pathBack.pop();
    path.push_back(currV);
      }
      // return vector
  return path;
}

// This function just displays the information of the person 1 and 2's building
// as well as the center building and all fo their nodes
void printStatementsIntial(BuildingInfo p1Build,
BuildingInfo p2Build, BuildingInfo closetBuilding,
long long node1, long long node2, long long nodec,
map<long long, Coordinates> Nodes) {
      cout << "Person 1's point:" << endl;
      cout << " "<<p1Build.Fullname << endl;
      cout << " (" << p1Build.Coords.Lat << ", " << p1Build.Coords.Lon << ")" << endl;
      cout << "Person 2's point:" << endl;
      cout << " " << p2Build.Fullname << endl;
      cout << " (" << p2Build.Coords.Lat << ", " << p2Build.Coords.Lon << ")"<< endl;
      cout << "Destination Building:" << endl;
      cout << " " << closetBuilding.Fullname << endl;
      cout << " (" << closetBuilding.Coords.Lat << ", " << closetBuilding.Coords.Lon << ")" << endl << endl;
      cout << "Nearest P1 node:" << endl;
      cout << node1 << endl;
      cout << " (" << Nodes[node1].Lat << ", " << Nodes[node1].Lon << ")"<<endl;
      cout << "Nearest P2 node:" << endl;
      cout << node2 << endl;
      cout << " (" << Nodes[node2].Lat << ", " << Nodes[node2].Lon << ")" << endl;
      cout << "Nearest destination node:" << endl;
      cout << nodec << endl;
      cout << " (" << Nodes[nodec].Lat << ", " << Nodes[nodec].Lon << ")" << endl;
}

// This function displays the path each person
// has to take to get to the center node
void printEnd(vector<long long> path1,
vector<long long> path2,
map<long long, double> distances1,
map<long long, double> distances2,
long long nodec, bool& pathFound) {
  cout << "Person 1's distance to dest: " << distances1[nodec] << " miles" << endl;
      cout << "Path: " << path1[0];
      for (int i = 1; i < path1.size(); i ++) {
        cout << "->" << path1[i];
      }
      cout << endl << endl;
      cout << "Person 2's distance to dest: " << distances2[nodec] << " miles" << endl;
      cout << "Path: " << path2[0];
      for (int i = 1; i < path2.size(); i ++) {
        cout << "->" << path2[i];
      }
      cout << endl << endl;
      // make pathfound true to get out of while loop
  pathFound = true;
}

// This is the function for when at least one person
// can't make it to the center building
// So we have to find a new building to go to
//There are a lot of things brought by reference here
void failureReprecussions(BuildingInfo closetBuilding,
BuildingInfo p1Build, BuildingInfo p2Build,
set<string>& excluded, vector<BuildingInfo> Buildings,
long long & nodec, vector<FootwayInfo> Footways,
map<long long, Coordinates> Nodes, map<long long, double> distances1,
map<long long,double> distances2) {
            cout << endl << "At least one person was unable to reach the destination building. Finding next closest building..." << endl << endl;
        // add the current closest building name in excluded set so we know not to visit there
        excluded.insert(closetBuilding.Fullname);
        // find a new closest building
        closetBuilding = nearestBuildingSet(p1Build, p2Build, Buildings, excluded);
        // find new closest node
        nodec = nearestNode(closetBuilding, Footways, Nodes);
        // display information
        cout << "New destination building:" << endl;
      cout << " " << closetBuilding.Fullname << endl;
      cout << " (" << closetBuilding.Coords.Lat << ", "
      << closetBuilding.Coords.Lon << ")" << endl << endl;
      cout << "Nearest destination node:" << endl;
      cout << nodec << endl;
      cout << " (" << Nodes[nodec].Lat << ", "
      << Nodes[nodec].Lon << ")" << endl;
      // check if this building still can't be reached
        if (distances1[nodec] >= INF || distances2[nodec] >= INF) {
          // if it can't be reached,
          // we then insert this buildings name in the excluded set
          excluded.insert(closetBuilding.Fullname);
        }
}

// This function checks whether or not either of the buildings exist
bool noNameChecker(BuildingInfo p1Build, BuildingInfo p2Build) {
        if (p1Build.Fullname == "") {
        cout << "Person 1's building not found" << endl;
        return true;
      } else if (p2Build.Fullname == "") {
        cout << "Person 2's building not found" << endl;
        return true;       }
  return false;
}

// This function is a call to the
// nearestNodes functions to create the nodes
void nodeCreater(long long& node1, long long& node2,
long long& nodec, BuildingInfo p1Build,
BuildingInfo p2Build, BuildingInfo closetBuilding,
vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes) {
       node1 = nearestNode(p1Build, Footways, Nodes);
       node2 = nearestNode(p2Build, Footways, Nodes);
       nodec = nearestNode(closetBuilding, Footways, Nodes);
}

// This is the implementation for the
// while loop in application ()
// and the calls to the different functions
void whileLoop(BuildingInfo p1Build, BuildingInfo p2Build,
map<long long, Coordinates> Nodes,
vector<FootwayInfo> Footways,
vector<BuildingInfo> Buildings,
graph<long long, double> G,
bool& pathFound) {
        // we check if the names of the buildings exist
        if (noNameChecker(p1Build, p2Build)) {
          // if not, then we exit the function
          // and we exit the while loop in applicaiton()
          pathFound = true;
        return;
      }
      // we find the closest building, and we send in an empty set
      // for the excluded set so that we don't deal with excluded names
      BuildingInfo closetBuilding = nearestBuildingSet(p1Build,
      p2Build, Buildings, {});
      // nodes for builings are created
      long long node1, node2, nodec;
      nodeCreater(node1, node2, nodec, p1Build, p2Build,
      closetBuilding, Footways, Nodes);
      // we print the intial statements for where
      // the nodes and buildings are located
      printStatementsIntial(p1Build, p2Build, closetBuilding,
      node1, node2, nodec, Nodes);
      // create a distances and predecessors map for each person's buidling
      // call dijkstra's to get shortest path
      map<long long, double> distances1;
      map<long long, long long> predecessors1;
      dijkstra(G, distances1, predecessors1, node1);
      map<long long, double> distances2;
      map<long long, long long> predecessors2;
      dijkstra(G, distances2, predecessors2, node2);
      // if we can't reach the other node from the origianal one
      // we just leave the loop
      if (distances1[node2] >= INF) {
          cout << "Sorry, destination unreachable." << endl << endl;
          pathFound = true;
          return;
        }
      set <string> excluded;
      // if we can't find a way to the current closest building, then we
      // find a new building in the following failureReprecussions call
      while (distances1[nodec] >= INF || distances2[nodec] >= INF) {
        failureReprecussions(closetBuilding,  p1Build,  p2Build, excluded,
        Buildings,  nodec, Footways, Nodes, distances1, distances2);
        }
        // once we find a suitable buidling we get the path
        vector <long long> path1 = getPath(predecessors1, nodec);
        vector <long long> path2 = getPath(predecessors2, nodec);
        // now we print the path we have found
      printEnd(path1,  path2,  distances1, distances2, nodec, pathFound);
}



/*
The creative component is to see how many miles a person would like to walk.
They start at a building a and then go to a building b.
If the distance they walked is less than the target distance, they are asked
where they want to go next (building c) to hit the target distance.
*/
void creative(map<long long, Coordinates>& Nodes,
vector<FootwayInfo>& Footways,
vector<BuildingInfo>& Buildings,
 graph<long long, double> G) {
    // we find how many miles they wanna walk, where they wanna start, and end
  double target = 0.0;
  cout << "How many miles would you like to walk today?" << endl;
  cin >> target;
  double currMiles = 0.0;
  string start;
  string end;
  cout << "Enter a starting location: " << endl;
  cin >> start;
  cout << "Enter an ending location: " << endl;
  cin >> end;
  // keep going until we hit our target distance
  while (currMiles < target) {
    // find start and ending buildings
  BuildingInfo startBuild = searchBuilding(start, Buildings);
  BuildingInfo endBuild = searchBuilding(end, Buildings);
  // while the starting and/or ending buildings dont exist,
  // prompt user for a new building
  while (startBuild.Fullname == "") {
    cout << "Starting building not found" << endl;
    cout << "Enter a new starting location: " << endl;
    cin >> start;
    startBuild = searchBuilding(start, Buildings);
  }
    while (endBuild.Fullname == "") {
    cout << "Ending building not found" << endl;
    cout << "Enter a new ending location: " << endl;
    cin >> end;
    endBuild = searchBuilding(end, Buildings);
  }
  // we find the nearest nodes to start and end buildings
  long long nodeS = nearestNode(startBuild, Footways, Nodes);
  long long nodeE = nearestNode(endBuild, Footways, Nodes);
  // make distance and predecessors map for start building
  map<long long, double> distances;
  map<long long, long long> predecessors;
  // dijkstra call
  dijkstra(G, distances, predecessors, nodeS);
  // if we can't reach end building form start,
  // then we tell user and prompt for new builing
    if (distances[nodeE] >= INF) {
      cout << "Sorry, no available path" << endl;
      cout << "Enter a new ending location: " << endl;
      cin >> end;
      // once we get a new building we go back to beginning of the loop
      continue;
    }
    // update how much we have walked
    currMiles += distances[nodeE];
    // if we haven't reach target, we inform user and repeat process
    if (currMiles < target) {
      cout << "You currently have " << currMiles << " miles" << endl;
      cout << "You still need " << target - currMiles
      << " miles to reach your target" << endl;
      cout << "Enter a new ending location: " << endl;
      string newEnd;
      cin >> newEnd;
      start = end;
      end = newEnd;
      // if we have reached target, then we leave loop
    } else {
      cout << "Congratulations, you have reached your target miles" << endl;
      break;
    }
  }
  return;
}

// This is the application for the myopenmaps,
// which will hold all the information
// and the calls for how to find the center building
// and the path between 2 others
void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double> G) {
  string person1Building, person2Building;
  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    // create the 2 person's building's from their names
    BuildingInfo p1Build = searchBuilding(person1Building, Buildings);
    BuildingInfo p2Build = searchBuilding(person2Building, Buildings);
     bool pathFound = false;
     // pathFound determine when to keep going in the loop
    while (!pathFound) {
      // call to whileLoop function to keep the application going
      whileLoop(p1Build,  p2Build,  Nodes, Footways, Buildings, G, pathFound);
      }
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
    }
  }





int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  graph <long long, double> G;
  // Milestone 5: add nodes to graph as vertexes
  for (auto n : Nodes) {
    G.addVertex(n.first);
  }
  // Milestone 6: add edges to graph-- follow what was said in the jumpstart
  for (auto fw : Footways) {
    for (int i = 0; i < fw.Nodes.size()-1; i++) {
      double lat1 = Nodes[fw.Nodes[i]].Lat;
      double long1 = Nodes[fw.Nodes[i]].Lon;
      double lat2 = Nodes[fw.Nodes[i+1]].Lat;
      double  long2 = Nodes[fw.Nodes[i+1]].Lon;
      double distance = distBetween2Points(lat1, long1, lat2, long2);
      G.addEdge(fw.Nodes[i], fw.Nodes[i+1], distance);
      G.addEdge(fw.Nodes[i+1], fw.Nodes[i], distance);
    }
  }

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


  //
  // TO DO: build the graph, output stats:
  //


  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative(Nodes,  Footways,  Buildings,  G);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
