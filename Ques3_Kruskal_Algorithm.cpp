#include <iostream>
#include <algorithm>
using namespace std;

struct Edge {
    int src, dest, weight;
};

int findParent(int node, int parent[]) {
    if (parent[node] == node)
        return node;
    return parent[node] = findParent(parent[node], parent);
}

void unionSets(int u, int v, int parent[], int rank[]) {
    int rootU = findParent(u, parent);
    int rootV = findParent(v, parent);

    if (rank[rootU] < rank[rootV]) {
        parent[rootU] = rootV;
    } else if (rank[rootU] > rank[rootV]) {
        parent[rootV] = rootU;
    } else {
        parent[rootV] = rootU;
        rank[rootU]++;
    }
}

bool compareEdges(const Edge& a, const Edge& b) {
    return a.weight < b.weight;
}

void kruskalMST(int V, int E, Edge edges[]) {
    sort(edges, edges + E, compareEdges);

    int parent[V], rank[V];
    for (int i = 0; i < V; i++) {
        parent[i] = i;
        rank[i] = 0;
    }

    Edge mst[V - 1];
    int mstIndex = 0;
    int totalWeight = 0;

    for (int i = 0; i < E; i++) {
        if (mstIndex == V - 1) break;

        int u = edges[i].src;
        int v = edges[i].dest;

        if (findParent(u, parent) != findParent(v, parent)) {
            mst[mstIndex++] = edges[i];
            totalWeight += edges[i].weight;
            unionSets(u, v, parent, rank);
        }
    }

    cout << "Edges in the Minimum Spanning Tree:\n";
    for (int i = 0; i < mstIndex; i++) {
        cout << mst[i].src << " -- " << mst[i].dest << " == " << mst[i].weight << endl;
    }
    cout << "Total weight of the MST: " << totalWeight << endl;
}

int main() {
    int V, E;
    cout << "Enter the number of vertices and edges: ";
    cin >> V >> E;

    Edge edges[E];
    cout << "Enter the edges (source, destination, weight):\n";
    for (int i = 0; i < E; i++) {
        cin >> edges[i].src >> edges[i].dest >> edges[i].weight;
    }

    kruskalMST(V, E, edges);

    return 0;
}

