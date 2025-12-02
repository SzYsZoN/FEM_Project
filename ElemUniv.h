#ifndef ELEMUNIV_H
#define ELEMUNIV_H

#include <vector>
#include <array>
#include <iostream>
#include "Gauss.h"

struct ElemUniv {
    std::vector<std::array<double,4>> dN_dXi;
    std::vector<std::array<double,4>> dN_dEta;

    explicit ElemUniv(int npc) {
        dN_dXi.resize(npc);
        dN_dEta.resize(npc);
    }

    void print() const {
        for (size_t p = 0; p < dN_dXi.size(); ++p) {
            std::cout << "GP " << p << "\n";
            std::cout << " dN/dXi:  ";
            for (int i = 0; i < 4; i++) std::cout << dN_dXi[p][i] << " ";
            std::cout << "\n dN/dEta: ";
            for (int i = 0; i < 4; i++) std::cout << dN_dEta[p][i] << " ";
            std::cout << "\n";
        }
    }
};

ElemUniv createElemUniv(const GaussQuadrature& quad);

#endif