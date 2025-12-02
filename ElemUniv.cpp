#include "ElemUniv.h"

ElemUniv createElemUniv(const GaussQuadrature& quad) {
    int npc = quad.points.size();
    ElemUniv eu(npc);

    for (int p = 0; p < npc; ++p) {
        double xi  = quad.points[p].xi;
        double eta = quad.points[p].eta;

        eu.dN_dXi[p][0] = -0.25 * (1 - eta);
        eu.dN_dXi[p][1] =  0.25 * (1 - eta);
        eu.dN_dXi[p][2] =  0.25 * (1 + eta);
        eu.dN_dXi[p][3] = -0.25 * (1 + eta);

        eu.dN_dEta[p][0] = -0.25 * (1 - xi);
        eu.dN_dEta[p][1] = -0.25 * (1 + xi);
        eu.dN_dEta[p][2] =  0.25 * (1 + xi);
        eu.dN_dEta[p][3] =  0.25 * (1 - xi);
    }

    return eu;
}
