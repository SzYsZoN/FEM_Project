#include <iostream>
#include <fstream>
#include "GlobalData.h"
#include "Grid.h"
#include "Gauss.h"
#include "ElemUniv.h"

static double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
    return s;
}

static std::vector<double> matvec(const std::vector<std::vector<double>>& A,
                                  const std::vector<double>& x)
{
    int n = (int)x.size();
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += A[i][j] * x[j];
        y[i] = s;
    }
    return y;
}

static std::vector<double> solveCG(const std::vector<std::vector<double>>& A,
                                   const std::vector<double>& b,
                                   double tol, int maxIter)
{
    int n = (int)b.size();
    std::vector<double> x(n, 0.0);  // start = 0

    std::vector<double> r = b;      // r = b - A*x = b
    std::vector<double> p = r;

    double rs_old = dot(r, r);
    if (std::sqrt(rs_old) < tol) return x;

    for (int it = 0; it < maxIter; it++)
    {
        std::vector<double> Ap = matvec(A, p);
        double denom = dot(p, Ap);
        if (std::fabs(denom) < 1e-30) break;

        double alpha = rs_old / denom;

        for (int i = 0; i < n; i++) x[i] += alpha * p[i];
        for (int i = 0; i < n; i++) r[i] -= alpha * Ap[i];

        double rs_new = dot(r, r);
        if (std::sqrt(rs_new) < tol) break;

        double beta = rs_new / rs_old;
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];

        rs_old = rs_new;
    }

    return x;
}


int main() {
    GlobalData data;
    data.npc = 2; // statyczne ustawienie liczby punktów całkowania 
    data.load("Test1.txt");
    data.print();


    GaussQuadrature gq(data.npc);       // 2x2 → 4 punkty
    ElemUniv eu = createElemUniv(gq);


    Grid grid(data.npc);
    grid.load("Test1.txt");
    
    std::ofstream fout("jacobians_output.txt");


    grid.assembleH(data,gq, eu);

    int n = data.nN;               
    double dt = data.SimulationStepTime;
    int steps = (int)(data.SimulationTime / dt);

    std::vector<double> T_old(n, data.InitialTemp);
    std::vector<double> T_new(n, 0.0);

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = grid.H_global[i][j] + grid.C_global[i][j] / dt;

    for (int s = 1; s <= steps; s++)
        {
        std::vector<double> b(n, 0.0);

        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += (grid.C_global[i][j] / dt) * T_old[j];

            b[i] = grid.P_global[i] + sum;
        }
        if (s == 1) {
        std::cout << "Iteration 0\n";
        // std::cout << "Matrix (H + C/dt)\n";
        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n; j++)
        //         std::cout << A[i][j] << " ";
        //     std::cout << "\n";
        // }

        // std::cout << "\nVector (P + (C/dt)*T0)\n";
        // for (int i = 0; i < n; i++)
        //     std::cout << b[i] << " ";
        // std::cout << "\n\n";
}

        // 1) rozwiąż A * T_new = b
        T_new = solveCG(A, b, 1e-10, 5000);   

        // 2) wypisz Tmin/Tmax
        double Tmin = T_new[0], Tmax = T_new[0];
        for (double v : T_new) { Tmin = std::min(Tmin, v); Tmax = std::max(Tmax, v); }

        std::cout << "t=" << s*dt << "  Tmin=" << Tmin << "  Tmax=" << Tmax << "\n";

        // 3) przejście do następnego kroku
        T_old.swap(T_new);
        }

    
    
    
    // przekierowanie cout → plik
    std::streambuf* oldCout = std::cout.rdbuf();
    std::cout.rdbuf(fout.rdbuf());

    for (auto& el : grid.elements){
        el.printAll(gq);
    }

    grid.printHGlobal();
    grid.printPGlobal();

    std::cout.rdbuf(oldCout); // przywróć cout do konsoli

    std::cout << "Zapisano wyniki do jacobians_output.txt\n";
    return 0;
}
