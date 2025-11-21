#include <iostream>
#include "GlobalData.h"
#include "Grid.h"

int main() {
    GlobalData data;
    data.load("Test3.txt");
    data.print();

    Grid grid;
    grid.load("Test3.txt");
    grid.print();

    return 0;
}
