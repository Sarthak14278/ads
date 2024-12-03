//Random_Quicksort

#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

int comparisons = 0; 

void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        comparisons++; 
        if (arr[j] <= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

int partition_r(int arr[], int low, int high) {
    srand(time(NULL));
    int random = low + rand() % (high - low + 1); 
    swap(&arr[random], &arr[high]);
    return partition(arr, low, high);
}

void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pi = partition_r(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++)
        cout << arr[i] << " ";
    cout << endl;
}

int main() {
    int n;
    cout << "Enter the size of the array: ";
    cin >> n;
    int arr[n];
    cout << "Enter the elements of the array: ";
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }

    quickSort(arr, 0, n - 1);

    cout << "Sorted array: ";
    printArray(arr, n);

    cout << "Number of comparisons: " << comparisons << endl;

    return 0;
}


//Select sort

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;


void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}


int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}


int partition_r(int arr[], int low, int high) {
    srand(time(NULL));
    int random = low + rand() % (high - low);
    swap(&arr[random], &arr[high]);
    return partition(arr, low, high);
}


int randomizedSelect(int arr[], int low, int high, int k) {
    if (low == high)
        return arr[low];

    int pi = partition_r(arr, low, high);
    int length = pi - low + 1;

    if (length == k)
        return arr[pi];
    else if (k < length)
        return randomizedSelect(arr, low, pi - 1, k);
    else
        return randomizedSelect(arr, pi + 1, high, k - length);
}

int main() {
    int n, k;
    cout << "Enter the size of the array: ";
    cin >> n;
    int arr[n];
    cout << "Enter the elements of the array: ";
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }
    cout << "Enter the value of k: ";
    cin >> k;

    int result = randomizedSelect(arr, 0, n - 1, k);
    cout << "The " << k << "-th smallest element is: " << result << endl;

    return 0;
}


//Kruskal algo

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

//Bellman Ford

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

//B Tree

#include <iostream>
using namespace std;

#define MIN_DEGREE 2

class BTreeNode {
public:
    int *keys;
    int t;
    BTreeNode **C;
    int n;
    bool leaf;

    BTreeNode(int _t, bool _leaf);
    void traverse();
    void insertNonFull(int key);
    void splitChild(int i, BTreeNode *y);
};

class BTree {
public:
    BTreeNode *root;
    int t;

    BTree(int _t) {
        root = NULL;
        t = _t;
    }

    void traverse() {
        if (root != NULL) root->traverse();
    }
    void insert(int key);
};

BTreeNode::BTreeNode(int _t, bool _leaf) {
    t = _t;
    leaf = _leaf;
    keys = new int[2 * t - 1];
    C = new BTreeNode *[2 * t];
    n = 0;
}

void BTreeNode::traverse() {
    for (int i = 0; i < n; i++) {
        if (!leaf) C[i]->traverse();
        cout << keys[i] << " ";
    }
    if (!leaf) C[n]->traverse();
}

void BTreeNode::insertNonFull(int key) {
    int i = n - 1;
    if (leaf) {
        while (i >= 0 && keys[i] > key) {
            keys[i + 1] = keys[i];
            i--;
        }
        keys[i + 1] = key;
        n++;
    } else {
        while (i >= 0 && keys[i] > key) i--;
        i++;
        if (C[i]->n == 2 * t - 1) {
            splitChild(i, C[i]);
            if (keys[i] < key) i++;
        }
        C[i]->insertNonFull(key);
    }
}

void BTreeNode::splitChild(int i, BTreeNode *y) {
    BTreeNode *z = new BTreeNode(y->t, y->leaf);
    z->n = t - 1;

    for (int j = 0; j < t - 1; j++)
        z->keys[j] = y->keys[j + t];

    if (!y->leaf) {
        for (int j = 0; j < t; j++)
            z->C[j] = y->C[j + t];
    }

    y->n = t - 1;

    for (int j = n; j >= i + 1; j--)
        C[j + 1] = C[j];

    C[i + 1] = z;

    for (int j = n - 1; j >= i; j--)
        keys[j + 1] = keys[j];

    keys[i] = y->keys[t - 1];
    n++;
}

void BTree::insert(int key) {
    if (root == NULL) {
        root = new BTreeNode(t, true);
        root->keys[0] = key;
        root->n = 1;
    } else {
        if (root->n == 2 * t - 1) {
            BTreeNode *s = new BTreeNode(t, false);
            s->C[0] = root;
            s->splitChild(0, root);
            int i = 0;
            if (s->keys[0] < key) i++;
            s->C[i]->insertNonFull(key);
            root = s;
        } else {
            root->insertNonFull(key);
        }
    }
}

int main() {
    BTree t(MIN_DEGREE);
    int choice, key;

    while (true) {
        cout << "\nB-Tree Menu\n";
        cout << "1. Insert key\n";
        cout << "2. Traverse B-tree\n";
        cout << "3. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "Enter key to insert: ";
                cin >> key;
                t.insert(key);
                break;
            case 2:
                cout << "Traversal of the B-tree is: ";
                t.traverse();
                cout << endl;
                break;
            case 3:
                cout << "Exiting...\n";
                return 0;
            default:
                cout << "Invalid choice. Please try again.\n";
        }
    }

    return 0;
}

//KMP

#include <iostream>
#include <vector>
using namespace std;

// Function to compute the Longest Prefix Suffix (LPS) array
void computeLPS(const string &pattern, vector<int> &lps) {
    int length = 0; // length of the previous longest prefix suffix
    lps[0] = 0;     // lps[0] is always 0
    int i = 1;

    while (i < pattern.length()) {
        if (pattern[i] == pattern[length]) {
            length++;
            lps[i] = length;
            i++;
        } else {
            if (length != 0) {
                // Backtrack the length
                length = lps[length - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

// Function to implement the KMP algorithm
void KMPAlgorithm(const string &text, const string &pattern) {
    int n = text.length();
    int m = pattern.length();

    // Create the LPS array
    vector<int> lps(m);
    computeLPSArray(pattern, lps);

    int i = 0; // index for text
    int j = 0; // index for pattern

    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }

        if (j == m) {
            cout << "Pattern found at index " << i - j << endl;
            j = lps[j - 1]; // Get the next position from LPS array
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1]; // Backtrack in the pattern
            } else {
                i++;
            }
        }
    }
}

int main() {
    string text, pattern;
    
    cout << "Enter the text: ";
    cin >> text;
    
    cout << "Enter the pattern to search: ";
    cin >> pattern;

    KMPAlgorithm(text, pattern);

    return 0;
}

//Tree
#include <iostream>
using namespace std;

// A B-Tree Node
class BTreeNode {
public:
    int* keys;          // Array of keys
    int t;              // Minimum degree
    BTreeNode** C;      // Array of child pointers
    int n;              // Current number of keys
    bool leaf;          // True if node is a leaf

    BTreeNode(int _t, bool _leaf);

    void traverse();

    BTreeNode* search(int k);

    void insertNonFull(int k);

    void splitChild(int i, BTreeNode* y);

    friend class BTree;
};

class BTree {
public:
    BTreeNode* root; // Pointer to root node
    int t;           // Minimum degree

    BTree(int _t) {
        root = nullptr;
        t = _t;
    }

    void traverse() {
        if (root != nullptr)
            root->traverse();
    }

    BTreeNode* search(int k) {
        return (root == nullptr) ? nullptr : root->search(k);
    }

    void insert(int k);
};

// Constructor for BTreeNode
BTreeNode::BTreeNode(int _t, bool _leaf) {
    t = _t;
    leaf = _leaf;
    keys = new int[2 * t - 1];
    C = new BTreeNode*[2 * t];
    n = 0;
}

// Traverse all nodes in a subtree rooted at this node
void BTreeNode::traverse() {
    int i;
    for (i = 0; i < n; i++) {
        if (!leaf)
            C[i]->traverse();
        cout << " " << keys[i];
    }
    if (!leaf)
        C[i]->traverse();
}

// Search a key in the subtree rooted at this node
BTreeNode* BTreeNode::search(int k) {
    int i = 0;
    while (i < n && k > keys[i])
        i++;

    if (keys[i] == k)
        return this;

    if (leaf)
        return nullptr;

    return C[i]->search(k);
}

// Insert a key into the B-Tree
void BTree::insert(int k) {
    if (root == nullptr) {
        root = new BTreeNode(t, true);
        root->keys[0] = k;
        root->n = 1;
    } else {
        if (root->n == 2 * t - 1) {
            BTreeNode* s = new BTreeNode(t, false);
            s->C[0] = root;
            s->splitChild(0, root);

            int i = 0;
            if (s->keys[0] < k)
                i++;
            s->C[i]->insertNonFull(k);

            root = s;
        } else {
            root->insertNonFull(k);
        }
    }
}

// Insert a key into a non-full node
void BTreeNode::insertNonFull(int k) {
    int i = n - 1;

    if (leaf) {
        while (i >= 0 && keys[i] > k) {
            keys[i + 1] = keys[i];
            i--;
        }

        keys[i + 1] = k;
        n++;
    } else {
        while (i >= 0 && keys[i] > k)
            i--;

        if (C[i + 1]->n == 2 * t - 1) {
            splitChild(i + 1, C[i + 1]);

            if (keys[i + 1] < k)
                i++;
        }
        C[i + 1]->insertNonFull(k);
    }
}

// Split the child y of this node
void BTreeNode::splitChild(int i, BTreeNode* y) {
    BTreeNode* z = new BTreeNode(y->t, y->leaf);
    z->n = t - 1;

    for (int j = 0; j < t - 1; j++)
        z->keys[j] = y->keys[j + t];

    if (!y->leaf) {
        for (int j = 0; j < t; j++)
            z->C[j] = y->C[j + t];
    }

    y->n = t - 1;

    for (int j = n; j >= i + 1; j--)
        C[j + 1] = C[j];

    C[i + 1] = z;

    for (int j = n - 1; j >= i; j--)
        keys[j + 1] = keys[j];

    keys[i] = y->keys[t - 1];
    n++;
}

// Main function to demonstrate the B-Tree
int main() {
    int t;
    cout << "Enter the minimum degree of the B-Tree: ";
    cin >> t;

    BTree tree(t);

    int choice, key;
    do {
        cout << "\n1. Insert\n2. Traverse\n3. Search\n4. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "Enter key to insert: ";
                cin >> key;
                tree.insert(key);
                break;
            case 2:
                cout << "B-Tree: ";
                tree.traverse();
                cout << endl;
                break;
            case 3:
                cout << "Enter key to search: ";
                cin >> key;
                if (tree.search(key) != nullptr)
                    cout << "Key " << key << " found in the tree.\n";
                else
                    cout << "Key " << key << " not found in the tree.\n";
                break;
            case 4:
                cout << "Exiting.\n";
                break;
            default:
                cout << "Invalid choice.\n";
        }
    } while (choice != 4);

    return 0;
}

//Suffix tree

#include <iostream>
#include <map>
#include <string>

class SuffixTreeNode {
public:
    std::map<char, SuffixTreeNode*> children; // Map to store child nodes
    int start;  // Start index of the edge
    int* end;   // End index of the edge (pointer for flexibility)

    SuffixTreeNode(int start, int* end) : start(start), end(end) {}
};

class SuffixTree {
private:
    SuffixTreeNode* root; // Root of the suffix tree
    std::string text;     // Input text
    int textLength;       // Length of the text

    // Helper function to print the suffix tree
    void printTree(SuffixTreeNode* node, int level = 0) {
        for (const auto& child : node->children) {
            std::cout << std::string(level * 4, ' '); // Indentation
            std::cout << "'" << text.substr(child.second->start, 
                        *(child.second->end) - child.second->start + 1) << "'\n";
            printTree(child.second, level + 1);
        }
    }

public:
    SuffixTree(const std::string& input) : text(input), textLength(input.length()) {
        root = new SuffixTreeNode(-1, new int(-1)); // Root node
        buildSuffixTree();
    }

    ~SuffixTree() {
        deleteTree(root);
    }

    // Function to delete the tree recursively
    void deleteTree(SuffixTreeNode* node) {
        for (auto& child : node->children) {
            deleteTree(child.second);
        }
        delete node->end;
        delete node;
    }

    // Function to build the suffix tree
    void buildSuffixTree() {
        for (int i = 0; i < textLength; ++i) {
            extendSuffixTree(i);
        }
    }

    // Extend the suffix tree for the ith suffix
    void extendSuffixTree(int pos) {
        SuffixTreeNode* currentNode = root;
        for (int i = pos; i < textLength; ++i) {
            char currentChar = text[i];
            if (currentNode->children.find(currentChar) == currentNode->children.end()) {
                // If the current character is not present, create a new node
                currentNode->children[currentChar] = new SuffixTreeNode(i, new int(textLength - 1));
            }
            currentNode = currentNode->children[currentChar];
        }
    }

    // Print the suffix tree
    void print() {
        std::cout << "Suffix Tree Structure:\n";
        printTree(root);
    }
};

int main() {
    std::string input = "banana";
    SuffixTree tree(input);
    tree.print();
    return 0;
}

