#include "GlobalData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
using namespace std;

static inline string trim(const string& s) { // usuwa biale znaki 
    size_t b = 0, e = s.size();
    while (b < e && isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && isspace(static_cast<unsigned char>(s[e-1]))) --e;
    return s.substr(b, e-b);
}
static inline string commas_to_spaces(string s) { //zamienia przecinki na spacje 
    for (char& c : s) if (c == ',') c = ' ';
    return s;
}

GlobalData::GlobalData()//konstruktor
    : SimulationTime(0), SimulationStepTime(0), Conductivity(0), Alfa(0),
      Tot(0), InitialTemp(0), Density(0), SpecificHeat(0), nN(0), nE(0),npc(0) {}//statyczne npc 

void GlobalData::load(const string& filename) {//wczytywanie danych z pliku txt 
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        return;
    }
    string line;
    while (getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '*') continue;
        string norm = commas_to_spaces(line);
        string k1, k2;
        stringstream ss(norm);
        if (!(ss >> k1)) continue;
        streampos pos = ss.tellg();
        if (ss >> k2) {
            string key2 = k1 + " " + k2;
            if (key2 == "Nodes number")     { ss >> nN; continue; }
            if (key2 == "Elements number")  { ss >> nE; continue; }
            ss.clear();
            ss.seekg(pos);
        }
        if      (k1 == "SimulationTime")          ss >> SimulationTime;
        else if (k1 == "SimulationStepTime")       ss >> SimulationStepTime;
        else if (k1 == "Conductivity")             ss >> Conductivity;
        else if (k1 == "Alfa")                     ss >> Alfa;
        else if (k1 == "Tot")                      ss >> Tot;
        else if (k1 == "InitialTemp")              ss >> InitialTemp;
        else if (k1 == "Density")                  ss >> Density;
        else if (k1 == "SpecificHeat")             ss >> SpecificHeat;
        else if (k1 == "nN")                       ss >> nN;
        else if (k1 == "nE")                       ss >> nE;
    }
}

void GlobalData::print() const {//wypisywanie danych globalnych ze struktury 
    cout << " GLOBAL DATA \n"
         << " SimulationTime     : " << SimulationTime     << "\n"
         << " SimulationStepTime : " << SimulationStepTime << "\n"
         << " Conductivity       : " << Conductivity       << "\n"
         << " Alfa               : " << Alfa               << "\n"
         << " Tot                : " << Tot                << "\n"
         << " InitialTemp        : " << InitialTemp        << "\n"
         << " Density            : " << Density            << "\n"
         << " SpecificHeat       : " << SpecificHeat       << "\n"
         << " Nodes number (hdr) : " << nN                 << "\n"
         << " Elements number(h) : " << nE                 << "\n"
         << " NPC (Gauss pts)    : " << npc                << "\n\n";
}
