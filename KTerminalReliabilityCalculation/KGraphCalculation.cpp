#include "KGraphCalculation.h"
#include "KGraphOperations.h"

int CountRec = 0;

double KFactoringProbabilityPST(Graph G, int variant)
{
    CountRec++;
    
    int e = 1;
    double p1 = 1;

    if (variant == 0) {
        if (KConnectivity(G) == false) return 0;
    } 
    if (TargetVertexQuality(G) == 1) return 1;
    else {
        //G = KParallelSeriesTransformation(G, p1);
        if (G.VertexCount == 2) {
            if (G.Targets[1] == 1 && G.Targets[2] == 1 && G.EdgeCount>0) return p1 * G.PArray[1];
            else return p1;
        }
        else if (G.VertexCount < 3) return p1;
    }
    
    int v = LastNotEmptyVertice(G);
    int u = G.FO[G.EdgeCount * 2 - 1];
    double p = G.PArray[G.EdgeCount * 2 - 1];
    //std::cout << CountRec << " ";
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
