#include "Gauss.h"
#include <vector>
#include <iostream>
#include <cmath>


    // Debug / print
void GaussQuadrature::print() const {
        std::cout << "Gauss points (" << points.size() << "):\n";
        for (size_t i = 0; i < points.size(); ++i) {
            std::cout << "  " << i+1 << ". "
                      << "xi=" << points[i].xi
                      << ", eta=" << points[i].eta
                      << ", w=" << points[i].weight
                      << "\n";
        }
    }
void GaussQuadrature::generate(int n) {
        switch (n) {
            case 1: generate1(); break;
            case 2: generate2(); break;
            case 3: generate3(); break;
            case 4: generate4(); break;

            default:
                std::cerr << "[GaussQuadrature] ERROR: Unsupported n = " << n << "\n";
        }
    }
void GaussQuadrature::generate1() {
        points = {
            {0.0, 0.0, 4.0}   // w = 2*2 bo 1D ma wagę 2
        };
    }

 void GaussQuadrature::generate2() {
        double a = 1.0 / std::sqrt(3.0);
        points = {
            { -a, -a, 1.0 },
            {  a, -a, 1.0 },
            {  a,  a, 1.0 },
            { -a,  a, 1.0 }
        };
    }

 void GaussQuadrature::generate3() {
        double a = std::sqrt(3.0/5.0);

        double w1 = 5.0 / 9.0;
        double w2 = 8.0 / 9.0;

        points = {
            // wiersz 1
            { -a, -a, w1*w1 }, { 0, -a, w2*w1 }, { a, -a, w1*w1 },
            // wiersz 2
            { -a,  0, w1*w2 }, { 0,  0, w2*w2 }, { a,  0, w1*w2 },
            // wiersz 3
            { -a,  a, w1*w1 }, { 0,  a, w2*w1 }, { a,  a, w1*w1 }
        };
    }

 void GaussQuadrature::generate4() {
        // standard Gauss-Legendre roots for 4-point
        double a = std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0));
        double b = std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0));

        double w1 = (18.0 - std::sqrt(30.0)) / 36.0;
        double w2 = (18.0 + std::sqrt(30.0)) / 36.0;

        double roots[4]   = { -a, -b,  b,  a };
        double weights[4] = { w1, w2, w2, w1 };

        points.clear();
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                points.emplace_back(roots[i], roots[j],
                                    weights[i] * weights[j]);
    }

