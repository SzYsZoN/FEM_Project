#include "ElemUniv.h"


ElemUniv createElemUniv(const GaussQuadrature& quad) {
    int npc = quad.points.size();
    ElemUniv eu(npc);
    

    for (int p = 0; p < npc; ++p) {
        double xi  = quad.points[p].xi;
        double eta = quad.points[p].eta;

        eu.N[p][0] = 0.25 * (1 - xi) * (1 - eta);
        eu.N[p][1] = 0.25 * (1 + xi) * (1 - eta);
        eu.N[p][2] = 0.25 * (1 + xi) * (1 + eta);
        eu.N[p][3] = 0.25 * (1 - xi) * (1 + eta);

        eu.dN_dXi[p][0] = -0.25 * (1 - eta);
        eu.dN_dXi[p][1] =  0.25 * (1 - eta);
        eu.dN_dXi[p][2] =  0.25 * (1 + eta);
        eu.dN_dXi[p][3] = -0.25 * (1 + eta);

        eu.dN_dEta[p][0] = -0.25 * (1 - xi);
        eu.dN_dEta[p][1] = -0.25 * (1 + xi);
        eu.dN_dEta[p][2] =  0.25 * (1 + xi);
        eu.dN_dEta[p][3] =  0.25 * (1 - xi);
    }

        int npc1d = std::sqrt(npc);   // dla 2x2 → npc = 4, npc1d = 2

    
    std::vector<double> ksi1D(npc1d);
    std::vector<double> w1D(npc1d);

    if (npc1d == 1) {
        ksi1D[0] = 0.0;
        w1D[0]   = 2.0;
    }
    else if (npc1d == 2) {
        double a = 1.0 / std::sqrt(3.0);
        ksi1D[0] = -a; ksi1D[1] = a;
        w1D[0] = 1.0;  w1D[1] = 1.0;
    }
    else if (npc1d == 3) {
        double a = std::sqrt(3.0/5.0);
        ksi1D[0] = -a; ksi1D[1] = 0.0; ksi1D[2] = a;
        w1D[0] = 5.0/9.0; w1D[1] = 8.0/9.0; w1D[2] = 5.0/9.0;
    }
    else if (npc1d == 4) {
        double a = std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0));
        double b = std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0));

        ksi1D[0] = -a; ksi1D[1] = -b; ksi1D[2] = b; ksi1D[3] = a;

        double w1 = (18.0 - std::sqrt(30.0)) / 36.0;
        double w2 = (18.0 + std::sqrt(30.0)) / 36.0;
        w1D[0] = w1; w1D[1] = w2; w1D[2] = w2; w1D[3] = w1;
    }


    for (int face = 0; face < 4; face++)
    {
        eu.surfaces[face].N.clear();
        eu.surfaces[face].detJ.clear();   
        eu.surfaces[face].w.clear();

        for (int p = 0; p < npc1d; p++)
        {
            double ksi, eta;

            switch (face)
            {
            case 0: // dół: eta = -1, ksi zmienne
                ksi = ksi1D[p];
                eta = -1.0;
                break;

            case 1: // prawa: ksi = +1, eta zmienne
                ksi =  1.0;
                eta = ksi1D[p];
                break;

            case 2: // góra: eta = +1, ksi zmienne
                ksi = ksi1D[p];
                eta =  1.0;
                break;

            case 3: // lewa: ksi = -1, eta zmienne
                ksi = -1.0;
                eta = ksi1D[p];
                break;
            }

            std::vector<double> N(4);
            N[0] = 0.25 * (1 - ksi) * (1 - eta);
            N[1] = 0.25 * (1 + ksi) * (1 - eta);
            N[2] = 0.25 * (1 + ksi) * (1 + eta);
            N[3] = 0.25 * (1 - ksi) * (1 + eta);

            eu.surfaces[face].N.push_back(N);
            eu.surfaces[face].w.push_back(w1D[p]);

        }
    }


    return eu;
}
