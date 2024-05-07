#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include "KGraphCalculation.h"
#include "KGraphOperations.h"
#include "TypesUnit.h"


Graph ReadingGraphFromFile();
void OutputInFile(double result);

int main()
{
    Graph G = ReadingGraphFromFile();
    double result = Probability(G);

    std::cout << result;
    OutputInFile(result);
}

Graph ReadingGraphFromFile() {
    Graph G;
    char comma;

    std::ifstream graphFile("Input.txt");

    if (graphFile.is_open()) {
        graphFile >> G.VertexCount;
        graphFile >> G.EdgeCount;

        for (int i{ 0 }; i < G.VertexCount + 1;i++) {
            graphFile >> G.KAO[i];
            if (i != G.VertexCount) graphFile >> comma;
        }
        for (int i{ 0 };i < G.EdgeCount * 2;i++) {
            graphFile >> G.FO[i];
            if (i != G.EdgeCount * 2 - 1) graphFile >> comma;
        }
        G.Targets[0] = 0;
        for (int i{ 1 }; i < G.VertexCount + 1;i++) {
            graphFile >> G.Targets[i];
            if (i != G.VertexCount) graphFile >> comma;
        }
        for (int i{ 0 }; i < 2 * G.EdgeCount; i++) G.PArray[i] = 0.5;
    }
    graphFile.close();
    return G;
}

void OutputInFile(double result) {
    std::ofstream out;
    out.open("Output.txt");
    if (out.is_open()) {
        out << result;
    }
    out.close();
    std::cout << "\nResult is saved to a file";
}