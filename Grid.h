#pragma once
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include "Jacobian.h"
#include "ElemUniv.h"

void printMatrix4(const std::vector<std::vector<double>>& M, const std::string& name);


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
    
    std::vector<std::array<double,4>> dNdx;
    std::vector<std::array<double,4>> dNdy;


    std::vector<std::vector<std::vector<double>>> H; // npc × 4 × 4
    std::vector<std::vector<double>> Hsum;           // 4 × 4

    Element(int id, Node* n1, Node* n2, Node* n3, Node* n4, int npc);
    void computeJacobian(const ElemUniv& eu);
    void computeH(const ElemUniv& eu, const GaussQuadrature& gq, double k);
    


    void print() const;
    void printH() const;
    void printAll(const GaussQuadrature& gq) const;
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
