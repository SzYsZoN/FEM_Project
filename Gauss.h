#ifndef GAUSS_H
#define GAUSS_H

#include <vector>

struct GaussPoint {
    double xi;
    double eta;
    double weight;

    GaussPoint() : xi(0), eta(0), weight(0) {}
    GaussPoint(double xi, double eta, double w)
        : xi(xi), eta(eta), weight(w) {}
};

class GaussQuadrature {
public:
    std::vector<GaussPoint> points;

    GaussQuadrature(int n){
        generate(n);
    }
    void print() const;

private:
    void generate(int n);
    void generate1();
    void generate2();
    void generate3();
    void generate4();
};

#endif
