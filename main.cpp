#include <iostream>
#include <fstream>
#include "GlobalData.h"
#include "Grid.h"
#include "Gauss.h"
#include "ElemUniv.h"


int main() {
    GlobalData data;
    data.npc = 1; // statyczne ustawienie liczby punktów całkowania 
    data.load("Test2.txt");
    data.print();


    GaussQuadrature gq(data.npc);       // 2x2 → 4 punkty
    ElemUniv eu = createElemUniv(gq);


    Grid grid(data.npc);
    grid.load("Test2.txt");
    
    std::ofstream fout("jacobians_output.txt");

    

    for (auto& el : grid.elements)
    el.computeJacobian(eu);

    // przekierowanie cout → plik
    std::streambuf* oldCout = std::cout.rdbuf();
    std::cout.rdbuf(fout.rdbuf());

    eu.print();
    grid.print();

    std::cout.rdbuf(oldCout); // przywróć cout do konsoli

    std::cout << "Zapisano wyniki do jacobians_output.txt\n";
    return 0;
}
