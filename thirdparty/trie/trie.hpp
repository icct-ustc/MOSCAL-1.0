#ifndef TRIE_H_
#define TRIE_H_

#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

class TrieNode {

public:

    // data
    TrieNode **child;
    int tier;
    int rank;

    // constructor
    TrieNode (): child(NULL), tier(-9527), rank(-9527) {}

    // constructor
    TrieNode (const int n_child): tier(-9527), rank(-9527) {
        if (n_child > 0) {
            child = new TrieNode*[n_child];
            for (int i=0; i<n_child; ++i) {
                child[i] = NULL;
            }
        } else {
            child = NULL;
        }
    }

    TrieNode (const int n_child, const int _tier, const int _rank): tier(_tier), rank(_rank) {
        if (n_child > 0) {
            child = new TrieNode*[n_child];
            for (int i=0; i<n_child; ++i) {
                child[i] = NULL;
            }
        } else {
            child = NULL;
        }
    }

    // destructor
   ~TrieNode () {
        delete [] child;
        child = NULL;
        tier = -9527;
        rank = -9527;
    }
};

class Trie {

public:

    // constructor
    Trie (): n_child(0), n_size(1) {
        root = new TrieNode();
    }

    // initialize
    void init (const int _n_child) {
        if (n_child==0 && n_size==1 && !root->child) {
            n_child = _n_child;
            root->child = new TrieNode*[n_child];
            for (int i=0; i<n_child; ++i) {
                root->child[i] = NULL;
            }
        } else {
            cout << "Error in initializing Trie!" << endl;
        }
    }

    // constructor
    Trie (const int _n_child): n_child(_n_child), n_size(1) {
        root = new TrieNode(n_child);
    }

    // constructor
    Trie (const Trie& rhs): n_child(rhs.n_child), n_size(rhs.n_size) {
        root = new TrieNode(n_child);
        for (int i=0; i<n_child; ++i) {
            if (rhs.root->child[i]) {
                root->child[i] = new TrieNode(n_child);
                clone(root->child[i],rhs.root->child[i]);
            } 
        }
    }

    // recursive clone
    void clone (TrieNode* dest, TrieNode* source) {
        dest->tier = source->tier;
        dest->rank = source->rank;
        for (int i=0; i<n_child-(source->tier); ++i) {
            if (source->child[i]) {
                dest->child[i] = new TrieNode(n_child-(source->tier)-i);
                clone(dest->child[i],source->child[i]);
            }
        }
    }

    // overload = operator
    Trie& operator= (const Trie& rhs) {
        if (this != &rhs) {
            n_child = rhs.n_child;
            n_size = rhs.n_size; 
            root = new TrieNode(n_child);
            for (int i=0; i<n_child; ++i) {
                if (rhs.root->child[i]) {
                    root->child[i] = new TrieNode(n_child);
                    clone(root->child[i],rhs.root->child[i]);
                } 
            }
        }
        return *this;
    }

   ~Trie () {
       for (int i=0; i<n_child; ++i) {
           if (root->child[i]) {
               clear(root->child[i]);
           }
       }
       delete root;
       root = NULL;
       n_child = n_size = 0;
    }

    // recursive clear
    void clear (TrieNode *p) {
       for (int i=0; i<n_child-(p->tier); ++i) {
           if (p->child[i]) {
               clear(p->child[i]);
           }
       }
       delete p;
       p = NULL;
    }

    // Try to insert a key-value pair; return true if succeed, false if exist;
    bool try_insert (const ivec& key, const int rank) {
        int cumtier = 0;
        TrieNode *p = root;
        for (unsigned int i=0; i<key.n_rows; ++i) {
            const int nk = key(i);
            cumtier += nk;
            if (!p->child[nk]) {
                p->child[nk] = new TrieNode(n_child-cumtier,cumtier,-9527);
                n_size += 1;
            }
            p = p->child[nk];
        }
        bool insertSucceed = false;
        if (p->rank<0) {
            insertSucceed = true;
            p->rank = rank;
        } 
        return insertSucceed;
    }

    // Return node pointer if key is in the trie.
    TrieNode* find (const ivec& key) const {
        TrieNode* p = root;
        for (unsigned int i=0; i<key.n_rows; ++i) {
            const int nk = key(i);
            if (!p->child[nk]) {
                return NULL;
            } else {
                p = p->child[nk];
            }
        }
        return p;
    }

    int size () const {
        return n_size;
    }

private:
    int n_child;
    int n_size;
    TrieNode* root;
};

#endif
