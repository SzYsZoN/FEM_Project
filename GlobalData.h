#pragma once
#include <string>


class GlobalData {
public:
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int nN;
    int nE;
    int npc;

    GlobalData();
    void load(const std::string& filename);
    void print() const;
};