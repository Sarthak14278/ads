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

