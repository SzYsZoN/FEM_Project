#include <iostream>
#include <fstream>
#include "GlobalData.h"
#include "Grid.h"
#include "Gauss.h"
#include "ElemUniv.h"



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
