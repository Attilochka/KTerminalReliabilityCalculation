#pragma once
const int NV = 25;
const int NE = 40;

struct Graph
{
    int KAO[NV + 1];
    int FO[NE * 2];
    int VertexCount;
    int EdgeCount;
    double PArray[NE * 2];
    int Targets[NV + 1];
};