#include <iostream>
#include <vector>
#include <limits.h>

using namespace std;

struct Edge {
    int source, destination, weight;
};

void bellmanFord(vector<Edge>& graph, int vertices, int source) {
    vector<int> distance(vertices, INT_MAX);
    distance[source] = 0;

    for (int i = 0; i < vertices - 1; ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            Edge e = graph[j];
            if (distance[e.source] != INT_MAX && distance[e.source] + e.weight < distance[e.destination]) {
                distance[e.destination] = distance[e.source] + e.weight;
            }
        }
    }

    for (size_t j = 0; j < graph.size(); ++j) {
        Edge e = graph[j];
        if (distance[e.source] != INT_MAX && distance[e.source] + e.weight < distance[e.destination]) {
            cout << "Graph contains negative weight cycle\n";
            return;
        }
    }

    cout << "Vertex Distance from Source\n";
    for (int i = 0; i < vertices; ++i) {
        cout << i << "\t" << distance[i] << "\n";
    }
}

int main() {
    int vertices, edges;
    cout << "Enter the number of vertices: ";
    cin >> vertices;
    cout << "Enter the number of edges: ";
    cin >> edges;

    vector<Edge> graph(edges);
    cout << "Enter each edge in the format (source destination weight):\n";
    for (int i = 0; i < edges; ++i) {
        cin >> graph[i].source >> graph[i].destination >> graph[i].weight;
    }

    int source;
    cout << "Enter the source vertex: ";
    cin >> source;

    bellmanFord(graph, vertices, source);

    return 0;
}

