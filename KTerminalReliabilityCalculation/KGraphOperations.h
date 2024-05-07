#pragma once
#include "TypesUnit.h"

bool KConnectivity(Graph G);
Graph DeleteEdgeGraph(Graph G, int u, int v);
Graph MergeVertex(Graph G, int u, int v);
int FindIndexEdgeForVertex(Graph G, int u, int v);
Graph KParallelSeriesTransformation(Graph G, double &p);
Graph Transformation(Graph G, int u, int v, int w, double &p);
bool KComponent(Graph &G);
Graph InducedKGraph(Graph G, int list[], int number);
