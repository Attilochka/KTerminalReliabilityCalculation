#pragma once
#include "TypesUnit.h"

bool KConnectivity(Graph G);
Graph DeleteEdgeGraph(Graph G, int u, int v);
Graph MergeVertex(Graph g, int u, int v);
int FindIndexEdgeForVertex(Graph G, int u, int v);
