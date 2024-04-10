#include "KGraphCalculation.h"
#include "KGraphOperations.h"

double KFactoringProbabilityPST(Graph G, int variant)
{
    int e = 1, p1 = 1;
    /* graph = DeleteEdgeGraph(graph, 1, 2);
     for (int i{ 0 }; i < graph.EdgeCount + 1; i++) {
         std::cout << graph.KAO[i] << " ";
     }
     graph = MergeVertex(graph, 1, 2);*/
    if (KConnectivity(G) == false && variant == 0) return 0;
    if (TargetVertexQuality(G) == 1) return 1;
    else {
        if (G.VertexCount == 3) {
            if (G.Targets[1] == 1 && G.Targets[2] == 1) return p1 * G.PArray[1];
            else return p1;
        }
        else if (G.VertexCount < 3) return p1;
    }
    int v = LastNotEmptyVertice(G);
    int u = G.FO[G.EdgeCount * 2 - 1];
    double p = G.PArray[G.EdgeCount * 2 - 1];


    return p1 * p * KFactoringProbabilityPST(MergeVertex(G, u, v), 1) + p1 * (1 - p) * KFactoringProbabilityPST(DeleteEdgeGraph(G, u, v), 0);
}

int TargetVertexQuality(Graph G)
{
    int result = 0;
    for (int i = 1; i < G.VertexCount + 1;i++) {
        if (G.Targets[i] == 1) {
            result++;
        }
    }
    return result;
}

unsigned int LastNotEmptyVertice(Graph G)
{
    int i = G.VertexCount;
    while (G.KAO[i] - G.KAO[i - 1] == 0) {
        i--;
    }
    return i;
}

double Probability(Graph G)
{
    double p = 1;
    return p * KFactoringProbabilityPST(G, 0);
}
