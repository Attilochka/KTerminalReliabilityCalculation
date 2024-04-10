#pragma once
const int NV = 20;
const int NE = 20;

struct Graph
{
    int KAO[NV + 1];
    int FO[NE * 2];
    int VertexCount;
    int EdgeCount;
    double PArray[NE * 2];
    int Targets[NV + 1];
};