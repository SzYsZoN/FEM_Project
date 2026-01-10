#include "Grid.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cctype>

using namespace std;

static inline string trim(const string& s) {
    size_t b = 0, e = s.size();
    while (b < e && isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && isspace(static_cast<unsigned char>(s[e-1]))) --e;
    return s.substr(b, e-b);
}
static inline string commas_to_spaces(string s) {
    for (char& c : s) if (c == ',') c = ' ';
    return s;
}

// ===== Node =====
Node::Node(int id, double x, double y) : id(id), x(x), y(y) {}
void Node::print() const {
    cout << "  - Node " << setw(3) << id << ": ("
         << fixed << setprecision(10) << x << ", " << y << ")\n";
}

//Wypisywanie macierzy 
void printMatrix4(const std::vector<std::vector<double>>& M, const std::string& name)
{
    std::cout << name << " =\n";
    for (int i = 0; i < 4; i++) {
        std::cout << "  ";
        for (int j = 0; j < 4; j++) {
            std::cout << std::setw(12) << M[i][j] << ' '; //wyrownanie do prawej 
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

// ===== Element =====
Element::Element(int id, Node* n1, Node* n2, Node* n3, Node* n4, int npc)
    : id(id), jac(npc)
{
    nodes[0] = n1;
    nodes[1] = n2;
    nodes[2] = n3;
    nodes[3] = n4;
}

void Element::computeJacobian(const ElemUniv& eu)
{
    int npc = eu.dN_dXi.size();


    jac.resize(npc);
    dNdx.resize(npc);
    dNdy.resize(npc);

    for (int p = 0; p < npc; p++) {

        // 1️⃣ Liczenie macierzy Jacobiego J
        double dx_dxi  = 0, dy_dxi = 0;
        double dx_deta = 0, dy_deta = 0;

        for (int i = 0; i < 4; i++) {
            dx_dxi  += nodes[i]->x * eu.dN_dXi[p][i];
            dy_dxi  += nodes[i]->y * eu.dN_dXi[p][i];

            dx_deta += nodes[i]->x * eu.dN_dEta[p][i];
            dy_deta += nodes[i]->y * eu.dN_dEta[p][i];
        }

        jac[p].J[0][0] = dx_dxi;
        jac[p].J[0][1] = dy_dxi;
        jac[p].J[1][0] = dx_deta;
        jac[p].J[1][1] = dy_deta;

        // 2️⃣ Liczenie detJ
        jac[p].detJ =
            jac[p].J[0][0] * jac[p].J[1][1] -
            jac[p].J[0][1] * jac[p].J[1][0];

        // 3️⃣ Liczenie odwrotności Jacobiego
        double invDet = 1.0 / jac[p].detJ;

        jac[p].invJ[0][0] =  jac[p].J[1][1] * invDet;
        jac[p].invJ[0][1] = -jac[p].J[0][1] * invDet;
        jac[p].invJ[1][0] = -jac[p].J[1][0] * invDet;
        jac[p].invJ[1][1] =  jac[p].J[0][0] * invDet;

        // 4️⃣ Liczenie pochodnych dN/dx, dN/dy
        for (int i = 0; i < 4; i++) {
            dNdx[p][i] =
                jac[p].invJ[0][0] * eu.dN_dXi[p][i] +
                jac[p].invJ[0][1] * eu.dN_dEta[p][i];

            dNdy[p][i] =
                jac[p].invJ[1][0] * eu.dN_dXi[p][i] +
                jac[p].invJ[1][1] * eu.dN_dEta[p][i];
        }
    }
}


void Element::computeH(const ElemUniv& eu, const GaussQuadrature& gq, double k)
{
    int npc = gq.points.size();

    H.resize(npc, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0)));
    Hsum.assign(4, std::vector<double>(4, 0.0));

    for (int p = 0; p < npc; p++)
    {
        double wx = gq.points[p].weight;  // dla 2x2 zawsze 1
        double detJ = jac[p].detJ;

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                double term =
                    dNdx[p][i] * dNdx[p][j] +
                    dNdy[p][i] * dNdy[p][j];

                H[p][i][j] = k * term * wx * detJ;

                Hsum[i][j] += H[p][i][j];
            }
        }
    }
}

void Element::printH() const
{
    std::cout << "==== Element " << id << " ====\n";

    // wypisz macierze dla każdego punktu całkowania
    for (int p = 0; p < H.size(); p++) {
        std::cout << "-- H for Gauss Point " << p << " --\n";
        printMatrix4(H[p], "H_pc");
    }

    // wypisz sumę wszystkich H_pc
    printMatrix4(Hsum, "Hsum (local stiffness matrix)");
}

void Element::computeHbc(const GlobalData& data, const ElemUniv& eu)
{
    Hbc.assign(4, std::vector<double>(4, 0.0));

    // Element nodes in global coordinates
    double x[4], y[4];
    for (int i = 0; i < 4; i++) {
        x[i] = nodes[i]->x;
        y[i] = nodes[i]->y;
    }

    // Pairs of local node indices for each face:
    // 0: bottom (N1-N2)
    // 1: right  (N2-N3)
    // 2: top    (N3-N4)
    // 3: left   (N4-N1)
    int faceNodes[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    for (int f = 0; f < 4; f++)
    {
        int a = faceNodes[f][0];
        int b = faceNodes[f][1];

        if (!nodes[a]->BC || !nodes[b]->BC) // warunek brzegowy jak nie to przechodzi do kolejnej ściany
            continue;

        
        double dx = x[b] - x[a];
        double dy = y[b] - y[a];
        double detJ = std::sqrt(dx*dx + dy*dy) / 2.0;

        for (int p = 0; p < eu.surfaces[f].N.size(); p++) 
        {
            const auto& N = eu.surfaces[f].N[p];
            double w = eu.surfaces[f].w[p];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Hbc[i][j] += data.Alfa * N[i] * N[j] * w * detJ;
                }
            }
        }
    }
}

void Element::computeP(const GlobalData& data, const ElemUniv& eu)
{
    P_local.assign(4, 0.0);

    double x[4], y[4];
    for (int i = 0; i < 4; i++) {
        x[i] = nodes[i]->x;
        y[i] = nodes[i]->y;
    }

    int faceNodes[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    for (int f = 0; f < 4; f++)
    {
        int a = faceNodes[f][0];
        int b = faceNodes[f][1];

        if (!nodes[a]->BC || !nodes[b]->BC) continue;

        double dx = x[b] - x[a];
        double dy = y[b] - y[a];
        double detJ = sqrt(dx*dx + dy*dy) / 2.0;

        for (int p = 0; p < eu.surfaces[f].N.size(); p++)
        {
            const auto& N = eu.surfaces[f].N[p];
            double w = eu.surfaces[f].w[p];

            for (int i = 0; i < 4; i++)
                P_local[i] += data.Alfa * data.Tot * N[i] * w * detJ;
        }
    }
}


void Element::print() const {
    cout << "--- Element " << id << " ---\n";

    for (size_t p = 0; p < jac.size(); p++) {
        cout << "Macierz Jakobiego dla punktu " << p << "\n";
        jac[p].print();

        cout << "dN/dx = ";
        for (double v : dNdx[p]) cout << v << ", ";
        cout << "\n";

        cout << "dN/dy = ";
        for (double v : dNdy[p]) cout << v << ", ";
        cout << "\n";
    }
}

void Element::printAll(const GaussQuadrature& gq) const
{
    std::cout << "\n=============================================\n";
    std::cout << "Obliczenia macierzy Jakobiego dla elementu nr " << id << "\n";

    int npc = gq.points.size();

    for (int p = 0; p < npc; p++)
    {
        std::cout << "Macierz Jakobiego dla " << p+1 << " punktu calkowania\n";

        std::cout << jac[p].J[0][0] << "  " << jac[p].J[0][1] << "\n"
                  << jac[p].J[1][0] << "  " << jac[p].J[1][1] << "\n";

        std::cout << "DetJ = " << jac[p].detJ 
                  << " dla " << p+1 << " punktu calkowania\n\n";

        // dNdx
        std::cout << "wartosc dN/dx rowna sie\n";
        for (int i = 0; i < 4; i++)
            std::cout << std::fixed << std::setprecision(6) << dNdx[p][i] << ", ";
        std::cout << "\n";

        // dNdy
        std::cout << "wartosc dN/dy rowna sie\n";
        for (int i = 0; i < 4; i++)
            std::cout << std::fixed << std::setprecision(6) << dNdy[p][i] << ", ";
        std::cout << "\n\n";
    }

    // ----------- Hsum (macierz lokalna) ----------
    std::cout << "H dla elementu - " << id << "\n";

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << std::setw(12) 
                      << std::fixed << std::setprecision(6)
                      << Hsum[i][j];
        }
        std::cout << "\n";
    }

    std::cout << "P dla elementu - " << id << "\n";

    for (int i = 0; i < 4; i++) {
        std::cout << std::setw(12) 
                  << std::fixed << std::setprecision(6)
                  << P_local[i];
        std::cout << "\n";
    }

    std::cout << "\n";
}


// ===== Grid =====
Grid::Grid(int npc) : nN(0), nE(0), npc(npc) {}

void Grid::load(const string& filename) {
    
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        return;
    }

    string line, section;
    while (getline(file, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (line.rfind("*Node", 0) == 0)    { section = "Node";    continue; }
        if (line.rfind("*Element", 0) == 0) { section = "Element"; continue; }
        if (line.rfind("*BC", 0) == 0)      { section = "BC";      continue; }

        if (section == "Node") {
            string norm = commas_to_spaces(line);
            stringstream ss(norm);
            int id; double x, y;
            if (ss >> id >> x >> y) nodes.emplace_back(id, x, y);
            continue;
        }

        if (section == "Element") {
            string norm = commas_to_spaces(line);
            stringstream ss(norm);
            int id, n1, n2, n3, n4;
            if (ss >> id >> n1 >> n2 >> n3 >> n4) {
                elements.emplace_back(
                    id,
                    &nodes[n1-1], &nodes[n2-1],
                    &nodes[n3-1], &nodes[n4-1],
                    npc
                );
            }
            continue;
        }

        if (section == "BC") {
            string norm = commas_to_spaces(line);
            stringstream ss(norm);
            int idx;
            while (ss >> idx){ 
                nodes[idx - 1].BC = true;
                bc.push_back(idx);
            }
            continue;
        }
    }

    nN = static_cast<int>(nodes.size());
    nE = static_cast<int>(elements.size());

    H_global.assign(nN, std::vector<double>(nN, 0.0)); // H

    P_global.assign(nN, 0.0); //P
}

void Grid::assembleH(const GlobalData& data,const GaussQuadrature& gq, const ElemUniv& eu)
{
    

    for (auto& el : elements) // liczenie lokalnych macierzy i agregacja
    {
        el.computeJacobian(eu);
        el.computeH(eu, gq, data.Conductivity);
        el.computeHbc(data, eu);
        el.computeP(data,eu);

        int ids[4] = {
            el.nodes[0]->id - 1,
            el.nodes[1]->id - 1,
            el.nodes[2]->id - 1,
            el.nodes[3]->id - 1
        };

        // H
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                H_global[ids[i]][ids[j]] += el.Hsum[i][j] + el.Hbc[i][j];
            }
        }
        // P 
        for (int i = 0; i < 4; i++)
        {    
            P_global[ids[i]] += el.P_local[i];
        }

        }
}


void Grid::print() const {
    cout << " GRID DATA \n";
    cout << "Nodes (actual)   : " << nN << "\n";
    cout << "Elements(actual) : " << nE << "\n";
    cout << "Nodes \n";
    for (const auto& node : nodes) node.print();

    cout << " Elements \n";
    for (const auto& el : elements) el.print();

    if (!bc.empty()) {
        cout << "  - BC nodes (" << bc.size() << "): ";
        for (size_t i = 0; i < bc.size(); ++i) {
            cout << bc[i] << (i + 1 < bc.size() ? ", " : "");
        }
        cout << "\n";
    }
}


void Grid::printHGlobal()const{
        for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            std::cout << H_global[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void Grid::printPGlobal() const {
    for(int i = 0; i < nN; i++)
        std::cout << P_global[i] << "\n";
}
