#include <iostream>

struct Jacobian {
    double J[2][2]{};
    double invJ[2][2]{};
    double detJ = 0.0;

    void print(int gp = -1) const {
        if (gp >= 0)
            std::cout << "Jacobian for Gauss Point " << gp << ":\n";
        else
            std::cout << "Jacobian:\n";

        std::cout << "  J = [ "
                  << J[0][0] << "  " << J[0][1] << "\n"
                  << "        "
                  << J[1][0] << "  " << J[1][1] << " ]\n";

        std::cout << "  detJ = " << detJ << "\n";

        std::cout << "  invJ = [ "
                  << invJ[0][0] << "  " << invJ[0][1] << "\n"
                  << "           "
                  << invJ[1][0] << "  " << invJ[1][1] << " ]\n";
    }
};
