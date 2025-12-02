#include "Grid.h"
#include <fstream>
#include <sstream>
#include <iomanip>
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
            while (ss >> idx) bc.push_back(idx);
            continue;
        }
    }

    nN = static_cast<int>(nodes.size());
    nE = static_cast<int>(elements.size());

   


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
