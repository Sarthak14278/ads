#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;

class SuffixTree {
public:
    struct Node {
        map<char, Node*> children;  // Children of the current node
        int start, *end;            // Start index and end index for the substring
        int suffixLink;             // Link to the deepest node with the same prefix
        Node(int start, int* end) : start(start), end(end), suffixLink(-1) {}
    };

    SuffixTree(const string& text) : text(text) {
        // Root of the tree
        root = new Node(-1, nullptr);
        build();
    }

    void build() {
        int n = text.size();
        int* end = new int(0);
        for (int i = 0; i < n; i++) {
            *end = i;
            extendSuffixTree(i, end);
        }
    }

    void print() {
        printTree(root, 0);
    }

private:
    string text;  // Input string
    Node* root;   // Root node of the suffix tree

    void extendSuffixTree(int pos, int* end) {
        // Here we would normally add the logic for extending the suffix tree.
        // In this naive implementation, we are just adding the suffixes one by one.
        cout << "Adding suffix: " << text.substr(pos) << endl;
    }

    void printTree(Node* node, int depth) {
        for (auto& child : node->children) {
            for (int i = 0; i < depth; i++) cout << "-";
            cout << child.first << " -> " << text.substr(child.second->start, *child.second->end - child.second->start + 1) << endl;
            printTree(child.second, depth + 1);
        }
    }
};

int main() {
    string text;
    cout << "Enter a string: ";
    cin >> text;

    SuffixTree tree(text);
    tree.print();  // Print the suffix tree (in a very simple way for now)

    return 0;
}
