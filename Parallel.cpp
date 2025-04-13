#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <omp.h>
#include<chrono>
#include<mpi.h>


#define TOTAL_THREADS 4
using namespace std;

using Graph = unordered_map<string, unordered_map<string, int>>;


struct Path {
    vector<string> nodes;
    int totalWeight;
    // Add a comparison operator
    bool operator>(const Path& other) const {
        return totalWeight > other.totalWeight;
    }
    bool operator==(const Path& other) const { 
        return (nodes == other.nodes); 
    }
};
struct Node {
    string id;
    int distance;
};

struct CompareNode {
    bool operator()(Node const& n1, Node const& n2) {
        // Nodes with smaller distances have higher priority
        return n1.distance > n2.distance;
    }
};

Path dijkstra(Graph& graph, string startNode, string endNode) {
    unordered_map<string, string> previousNodes;
    unordered_map<string, int> distances;
    priority_queue<Node, vector<Node>, CompareNode> nodesQueue;


    // Initialize distances to infinity and previous nodes to empty
	    for (auto& pair : graph) {
		distances[pair.first] = INT_MAX;
		previousNodes[pair.first] = "";
	    }

    // The distance from the start node to itself is 0
    distances[startNode] = 0;
    nodesQueue.push({ startNode, 0 });

    while (!nodesQueue.empty()) {
        Node currentNode = nodesQueue.top();
        nodesQueue.pop();

        // For each neighbor of the current node
        
        for (auto& pair : graph[currentNode.id]) {
            string neighbor = pair.first;
            int weight = pair.second;

            // If a shorter path to the neighbor has been found
            if (distances[currentNode.id] + weight < distances[neighbor]) {
                distances[neighbor] = distances[currentNode.id] + weight;
                previousNodes[neighbor] = currentNode.id;
                nodesQueue.push({ neighbor, distances[neighbor] });
            }
        }
    }

    // Build the shortest path from startNode to endNode
    Path shortestPath;
    shortestPath.totalWeight = distances[endNode];
    for (string node = endNode; node != ""; node = previousNodes[node]) {
        shortestPath.nodes.push_back(node);
    }
    reverse(shortestPath.nodes.begin(), shortestPath.nodes.end());

    return shortestPath;
}

// ... rest of the code ...

vector<Path> yenKSP(Graph& graph, string startNode, string endNode, int K) {
    vector<Path> kShortestPaths;
    vector<Path> ShortestPaths;
    priority_queue<Path, vector<Path>, greater<Path>> candidatePaths;

    // Initialize with the shortest path from startNode to endNode
    Path shortestPath = dijkstra(graph, startNode, endNode);
    candidatePaths.push(shortestPath);

    while (!candidatePaths.empty() && kShortestPaths.size() < K) {
        Path currPath = candidatePaths.top();
        candidatePaths.pop();

        // Add the current path to the shortest paths list
        kShortestPaths.push_back(currPath);

        // Parallelize the loop over each node in currPath
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < currPath.nodes.size() - 1; ++i) {
            string spurNode = currPath.nodes[i];
            Graph graphCopy;
            Path spurPath;
            Path newPath;

            // Critical section to protect the graph copy
            #pragma omp critical
            {
                graphCopy = graph;
                // Remove the nodes and edges of the current path and the previously found shortest paths
                for (size_t j = 0; j <= i; ++j) {
                    string node = currPath.nodes[j];
                    if (node != spurNode) {
                        graphCopy.erase(node);
                    }
                }
                for (const Path& path : kShortestPaths) {
                    if (i < path.nodes.size() && path.nodes[i] == spurNode) {
                        graphCopy[spurNode].erase(path.nodes[i + 1]);
                    }
                }
            }

            // Calculate the shortest path from spurNode to endNode
            spurPath = dijkstra(graphCopy, spurNode, endNode);

             // Create a new path that combines the current path and the spur path
		 //Path newPath;
		 //Path newPath1;
		 bool a = true;
		 int index = 0;
		 for (int j = 0; j <= i; j++)
		 {
		     if (currPath.nodes[j] == spurNode)
		     {
			 a = false;
			 index = j;
			 break;
		     }
		 }
		 if(a && index==0)
		 newPath.nodes.insert(newPath.nodes.end(), currPath.nodes.begin(), currPath.nodes.begin() + i + 1);
		 if (!a && index != 0)
		     if (currPath.nodes.size() > 1) {
			 //newPath.nodes.push_back(currPath.nodes[0]);
			 for (int c = 0; c < index; c++)
			 {
			    newPath.nodes.push_back(currPath.nodes[c]);
			 }
		     }
		 newPath.nodes.insert(newPath.nodes.end(), spurPath.nodes.begin(), spurPath.nodes.end());

		 /*if (currPath.nodes[0] != spurNode)
		 newPath.totalWeight = currPath.totalWeight - graph[currPath.nodes[i]][currPath.nodes[i + 1]] + spurPath.totalWeight;
		 else*/
		 newPath.totalWeight = spurPath.totalWeight;

		 // Add the new path to candidate paths
		 newPath.totalWeight += index;
            // Critical section to protect the shared data structures
            #pragma omp critical
            {
                // Check if newPath is already in ShortestPaths
                auto it = std::find(ShortestPaths.begin(), ShortestPaths.end(), newPath);
                if (it == ShortestPaths.end() && newPath.totalWeight > 0) {
                    ShortestPaths.push_back(newPath);
                    candidatePaths.push(newPath);
                }
            }
        }
    }

    return kShortestPaths;
}


// ... rest of the code ...






void printPaths(const vector<Path>& paths) {
    int pathss = 0;
    for (int i = 0; i < paths.size(); ++i) {
        if (paths[i].totalWeight > 0)
        {
            cout << "Path " << pathss + 1 << ": ";
            for (const string& node : paths[i].nodes) {
                cout << node << " -> ";
            }
            cout << "END, Total weight: " << paths[i].totalWeight << endl;
            pathss++;
        }
    }
}

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

    // Get the rank of the current process and the total number of processes
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    omp_set_num_threads(TOTAL_THREADS);
    
    ifstream file("/home/azeem/Downloads/doctorwho.csv");
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file." << endl;
        return 1;
    }

    Graph graph;
    unordered_set<string> uniqueNodes;
    

    string line;
    getline(file, line); // Skip header line
    
	while (getline(file, line)) {
	stringstream ss(line);
	string source, target, weightStr;
	getline(ss, source, ',');
	getline(ss, target, ',');
	getline(ss, weightStr, ',');
	int weight =0;
	try {
	    weight = stoi(weightStr);
	} catch (const std::invalid_argument& e) {
	    cerr << "Invalid argument: " << weightStr << endl;
	}

	// Add edge from source to target with weight
	graph[source][target] = weight;

	// Check if the source and target nodes are new
	uniqueNodes.insert(source);
	uniqueNodes.insert(target);
	}
    



    /*string line;
    getline(file, line); // Skip header line
    while (getline(file, line)) {
        stringstream ss(line);
        string source, target, weightStr;
        getline(ss, source, ',');
        getline(ss, target, ',');
        getline(ss, weightStr, ',');
        int weight = stoi(weightStr);

        // Add edge from source to target with weight
        graph[source][target] = weight;

        // Check if the source and target nodes are new
        uniqueNodes.insert(source);
        uniqueNodes.insert(target);
    }*/

    file.close();

    //Number of unique nodes in the graph
    int N = uniqueNodes.size();

        /*ifstream file("C:\\Users\\Ahsan\\Downloads\\email-EuAll.txt\\Email-EuAll.txt");
        if (!file.is_open()) {
            cerr << "Error: Unable to open the file." << endl;
            return 1;
        }

        Graph graph;
        unordered_set<string> uniqueNodes;

        string line;
        // Skip the header lines
        for (int i = 0; i < 4; ++i) {
            getline(file, line);
        }

        while (getline(file, line)) {
            stringstream ss(line);
            string source, target;
            getline(ss, source, '\t');
            getline(ss, target);

            // Assign a default weight of 1 to every edge
            int weight = 1;

            // Add edge from source to target with weight
            graph[source][target] = weight;

            // Check if the source and target nodes are new
            uniqueNodes.insert(source);
            uniqueNodes.insert(target);
        }

        file.close();

        // Number of unique nodes in the graph
        int N = uniqueNodes.size();*/

    // Calculate the k shortest paths using Yen's algorithm
    // Start measuring time
    
    
    cout<<"Parallel Version"<<endl;
    
    for(int i=1;i<=10;i++){
    if(i%size==rank){
    srand(time(NULL));
    int random_number1 = rand() % N + 1;  // N is the number of unique nodes
    int random_number2 = rand() % N + 1;
    auto it1 = uniqueNodes.begin();
	    auto it2 = uniqueNodes.begin();
	    std::advance(it1, random_number1 - 1);
	    std::advance(it2, random_number2 - 1);
	    string startNode = *it1;
	    string endNode = *it2;
	    
	    cout<<endl;
	    cout<<"Execution Number: "<<i+1<<endl;
	    cout<<"Start Node: "<<startNode<<endl;
	    cout<<"End Node: "<<endNode<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    vector<Path> kShortestPaths = yenKSP(graph, startNode, endNode, 2);

    // Stop measuring time and calculate the elapsed time
    auto finish = std::chrono::high_resolution_clock::now();

    // Print the k shortest paths
    printPaths(kShortestPaths);
    
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	cout<<endl;
	}
	
	}
	MPI_Finalize();
    return 0;
}