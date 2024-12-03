#include <iostream>
using namespace std;

// Class to represent a node in the binary search tree
class TreeNode {
public:
    int data;
    TreeNode* left;
    TreeNode* right;

    TreeNode(int value) : data(value), left(NULL), right(NULL) {}
};

// Class to implement the binary search tree
class BinarySearchTree {
private:
    TreeNode* root;

    // Helper function to insert a value into the tree
    TreeNode* insertNode(TreeNode* node, int value) {
        if (node == NULL) {
            return new TreeNode(value);
        }
        if (value < node->data) {
            node->left = insertNode(node->left, value);
        } else if (value > node->data) {
            node->right = insertNode(node->right, value);
        }
        return node;
    }

    // Helper function to search for a value in the tree
    bool searchNode(TreeNode* node, int value) {
        if (node == NULL) {
            return false;
        }
        if (node->data == value) {
            return true;
        }
        if (value < node->data) {
            return searchNode(node->left, value);
        } else {
            return searchNode(node->right, value);
        }
    }

    // Helper function for in-order traversal
    void inOrder(TreeNode* node) {
        if (node != NULL) {
            inOrder(node->left);
            cout << node->data << " ";
            inOrder(node->right);
        }
    }

public:
    // Constructor
    BinarySearchTree() : root(NULL) {}

    // Method to insert a value into the tree
    void insert(int value) {
        root = insertNode(root, value);
    }

    // Method to search for a value in the tree
    bool search(int value) {
        return searchNode(root, value);
    }

    // Method to display the tree in sorted order
    void display() {
        inOrder(root);
        cout << endl;
    }
};

// Main function
int main() {
    BinarySearchTree bst;

    cout << "Binary Search Tree Operations:\n";
    while (true) {
        cout << "\n1. Insert\n2. Search\n3. Display\n4. Exit\nEnter your choice: ";
        int choice, value;
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "Enter value to insert: ";
                cin >> value;
                bst.insert(value);
                cout << "Value inserted.\n";
                break;
            case 2:
                cout << "Enter value to search: ";
                cin >> value;
                if (bst.search(value)) {
                    cout << "Value found in the tree.\n";
                } else {
                    cout << "Value not found in the tree.\n";
                }
                break;
            case 3:
                cout << "Tree elements in sorted order: ";
                bst.display();
                break;
            case 4:
                cout << "Exiting program.\n";
                return 0;
            default:
                cout << "Invalid choice. Please try again.\n";
        }
    }
}

// insertion -> best O(log N) and worst O(N)
// search -> same to same 
// inorder traveral -> O(N)
// space comp -> samme 
// preformance degrade to o of n when unbalanced

