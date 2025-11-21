#pragma once
#include <vector>
#include <string>
#include <iostream>

class Node {
public:
    int id;
    double x, y;
    Node(int id, double x, double y);
    void print() const;
};

class Element {
public:
    int id;
    Node* nodes[4];          // wskaźniki do węzłów
    Element(int id, Node* n1, Node* n2, Node* n3, Node* n4);
    void print() const;
};

class Grid {
public:
    int nN;
    int nE;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> bc;

    Grid();
    void load(const std::string& filename);
    void print() const;
};
