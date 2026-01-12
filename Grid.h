#pragma once
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include "Jacobian.h"
#include "ElemUniv.h"
#include "GlobalData.h"

void printMatrix4(const std::vector<std::vector<double>>& M, const std::string& name);


class Node {
public:
    int id;
    double x, y;
    bool BC = false;
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


    std::vector<std::vector<std::vector<double>>> H; 
    std::vector<std::vector<double>> Hsum;           

    std::vector<std::vector<double>> Hbc; 

    std::vector<double> P_local; 

    std::vector<std::vector<std::vector<double>>> C; 
    std::vector<std::vector<double>> Csum;           


    Element(int id, Node* n1, Node* n2, Node* n3, Node* n4, int npc);
    void computeJacobian(const ElemUniv& eu);
    void computeH(const ElemUniv& eu, const GaussQuadrature& gq, double k);
    void computeHbc(const GlobalData& data, const ElemUniv& eu);
    void computeP(const GlobalData& data, const ElemUniv& eu);
    void computeC(const ElemUniv& eu, const GaussQuadrature& gq, const GlobalData& data);
    


    void print() const;
    void printH() const;
    void printAll(const GaussQuadrature& gq) const;
};

class Grid {
public:

    std::vector<std::vector<double>> H_global;
    std::vector<double> P_global;
    std::vector<std::vector<double>> C_global;

    int npc;
    int nN;
    int nE;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> bc;

    Grid(int npc);
    void load(const std::string& filename);
    void assemble(const GlobalData& data,const GaussQuadrature& gq, const ElemUniv& eu);
    void print() const;
    void printHGlobal()const;
    void printPGlobal() const;
    void printCGlobal() const;
};
