
#include <vector>
#include <string>
#include <iostream>
#include "Jacobian.h"
#include "ElemUniv.h"

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
    Node* nodes[4];
    std::vector<Jacobian> jac;

    Element(int id, Node* n1, Node* n2, Node* n3, Node* n4, int npc);
    void computeJacobian(const ElemUniv& eu);
    std::vector<std::array<double,4>> dNdx;
    std::vector<std::array<double,4>> dNdy;


    void print() const;
};

class Grid {
public:
    int npc;
    int nN;
    int nE;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> bc;

    Grid(int npc);
    void load(const std::string& filename);
    void print() const;
};
