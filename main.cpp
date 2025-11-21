#include <iostream>
#include "GlobalData.h"

int main() {
    GlobalData data;
    data.load("Test1.txt");
    data.print();
    return 0;
}
