#include "KGraphOperations.h"
#include "KGraphCalculation.h"

bool KConnectivity(Graph G)
{
    int sum = 1, sum1 = 0, l = 1, LB;
    bool Result;
    int Spot[NV + 1], B[NV + 1], A[NV + 1];
    //std::vector<int> B(0);

    for (int i{ 1 }; i < G.VertexCount + 1; i++) {
        if (G.Targets[i] == 1) {
            A[0] = i;
            break;
        }
    }
    for (int i{ 1 }; i < G.VertexCount + 1; i++) Spot[i] = 0;
    Spot[A[0]] = 1;
    while (l > 0) {
        LB = 0;
        //B.resize(0);
        for (int i{ 0 };i < l;i++) {
            for (int j{ G.KAO[A[i] - 1] }; j < G.KAO[A[i]]; j++) {
                if (Spot[G.FO[j]] == 0) {
                    //B.resize(B.size()+1);
                    /* std::cout << B.size()*/;
                    //B[B.size() - 1] = G.FO[j];];
                    LB++;
                    B[LB - 1] = G.FO[j];
                    Spot[G.FO[j]] = 1;
                    if (G.Targets[G.FO[j]] == 1) sum++;
                }

            }
        }
        /*A.resize(B.size());
        l = A.size();*/
        l = LB;
        if (l > 0) {
            for (int i{ 0 }; i < l;i++) {
                A[i] = B[i];
            }
        }
    }


    sum1 = 0;
    for (int i{ 1 }; i < G.VertexCount + 1;i++) {
        if (G.Targets[i] == 1) sum1++;
        if (sum == sum1) {
            Result = true;
        }
        else Result = false;
        if (G.VertexCount == 1) Result = true;
    }
    return Result;
}

Graph DeleteEdgeGraph(Graph G, int u, int v)
{
    Graph newGraph;
    int indexU = 0, indexV;
    indexU = FindIndexEdgeForVertex(G, u, v);
    indexV = FindIndexEdgeForVertex(G, v, u);
    //std::cout << indexU << " " << indexV<< "\n";
    if (indexU < G.EdgeCount * 2 && indexV < G.EdgeCount * 2) {
        newGraph.VertexCount = G.VertexCount;
        newGraph.EdgeCount = G.EdgeCount - 1;
        if (indexU < indexV) {
            for (int i{ 0 }; i < u; i++) {
                newGraph.KAO[i] = G.KAO[i];
            }
            for (int i{ u }; i < v; i++) {
                newGraph.KAO[i] = G.KAO[i] - 1;
            }
            for (int i{ v }; i <= G.VertexCount; i++) {
                newGraph.KAO[i] = G.KAO[i] - 2;
            }
            if (indexU > 0) {
                for (int i{ 0 }; i < indexU; i++) {
                    newGraph.FO[i] = G.FO[i];
                    newGraph.PArray[i] = G.PArray[i];
                }
            }
            if (indexV > indexU + 1) {
                for (int i{ indexU }; i < indexV - 1; i++) {
                    newGraph.FO[i] = G.FO[i + 1];
                    newGraph.PArray[i] = G.PArray[i + 1];
                }
            }
            if (G.EdgeCount * 2 > indexV + 1) {
                for (int i{ indexV - 1 }; i < G.EdgeCount * 2 - 2; i++) {
                    newGraph.FO[i] = G.FO[i + 2];
                    newGraph.PArray[i] = G.PArray[i + 2];
                }
            }
        }
        else {
            for (int i{ 0 }; i < v - 1; i++) {
                newGraph.KAO[i] = G.KAO[i];
            }
            for (int i{ v }; i < u - 1; i++) {
                newGraph.KAO[i] = G.KAO[i];
            }
            for (int i{ u }; i <= G.VertexCount; i++) {
                newGraph.KAO[i] = G.KAO[i] - 2;
            }
            if (indexV > 0) {
                for (int i{ 0 }; i < indexV - 2; i++) {
                    newGraph.FO[i] = G.FO[i];
                    newGraph.PArray[i] = G.PArray[i];
                }
            }
            if (indexU > indexV + 1) {
                for (int i{ indexV }; i < indexU - 1;i++) {
                    newGraph.FO[i] = G.FO[i + 1];
                    newGraph.PArray[i] = G.PArray[i + 1];
                }
            }
            if (G.EdgeCount * 2 > indexU + 1) {
                for (int i{ indexU - 1 }; i < G.EdgeCount * 2 - 3; i++) {
                    newGraph.FO[i] = G.FO[i + 2];
                    newGraph.PArray[i] = G.PArray[i + 2];
                }
            }
        }
    }
    for (int i{ 0 };i < newGraph.VertexCount + 1; i++) {
        newGraph.Targets[i] = G.Targets[i];
    }
    return newGraph;
}

Graph MergeVertex(Graph G, int u, int v)
{
    int KAO[NV + 1], FO[2 * NE], S[NV + 1];
    double PArray[2 * NE];
    Graph Result;
    int i, j, l = 0, k, boolka;

    for (i = 1; i < v; i++) {
        S[i] = i;
    }
    S[v] = u;

    for (i = v + 1; i < G.VertexCount + 1; i++) S[i] = i - 1;

    KAO[0] = 0;

    for (i = 1; i < G.VertexCount + 1; i++)
    {
        if (i == u) { boolka = 1; }
        else
        {
            if (i == v) { boolka = 2; }
            else
            {
                k = 0;
                for (j = G.KAO[i - 1]; j < G.KAO[i]; j++)
                {
                    if ((G.FO[j] == u) || (G.FO[j] == v)) { k++; }
                }
                if (k == 2) { boolka = 3; }
                else { boolka = 4; }
            }
        }

        switch (boolka)
        {
        case 1:
            KAO[S[i]] = KAO[S[i] - 1];
            for (j = G.KAO[i - 1]; j < G.KAO[i]; j++)
            {
                if (G.FO[j] != v)
                {
                    ++KAO[S[i]];
                    FO[l] = S[G.FO[j]];
                    // k=SearchEdge(G,cutv2,G.FO[j]);
                    bool b = true;
                    for (int m = G.KAO[v - 1]; m < G.KAO[v]; m++)
                    {
                        if (G.FO[m] == G.FO[j]) {
                            b = false;
                            k = m;
                        }
                    }if (b == true) { k = G.EdgeCount * 2; }
                    //
                    if (k != G.EdgeCount * 2)
                    {
                        PArray[l] = G.PArray[j] + G.PArray[k] - G.PArray[j] * G.PArray[k];
                    }
                    else { PArray[l] = G.PArray[j]; }
                    ++l;
                }
            }
            for (j = G.KAO[v - 1]; j < G.KAO[v]; j++)
            {
                if (G.FO[j] != u)
                {//if (SearchEdge(G,cutv1,G.FO[j])==G.EdgNumb*2)
                    bool b = true;
                    int o;
                    for (int m = G.KAO[u - 1]; m < G.KAO[u]; m++)
                    {
                        if (G.FO[m] == G.FO[j]) {
                            b = false;
                            o = m;
                        }
                    }if (b == true) { o = G.EdgeCount * 2; } if (o == G.EdgeCount * 2)
                        //
                    {
                        ++KAO[S[i]];
                        FO[l] = S[G.FO[j]];
                        PArray[l] = G.PArray[j];
                        ++l;
                    }
                }
            }
            break;
        case 2:
            break;
        case 3:

            KAO[S[i]] = KAO[S[i] - 1] + G.KAO[i] - G.KAO[i - 1] - 1;
            for (j = G.KAO[i - 1]; j < G.KAO[i]; j++)
            {
                if (G.FO[j] != v)
                {
                    if (G.FO[j] != u)
                    {
                        FO[l] = S[G.FO[j]];
                        PArray[l] = G.PArray[j];
                        ++l;
                    }
                    else
                    {
                        FO[l] = S[G.FO[j]];
                        for (k = G.KAO[i - 1]; k < G.KAO[i]; k++)
                        {
                            if (G.FO[k] == v)
                            {
                                PArray[l] = G.PArray[j] + G.PArray[k] - G.PArray[j] * G.PArray[k];
                            }
                        }
                        ++l;
                    }
                }
            }
            break;
        case 4:

            KAO[S[i]] = KAO[S[i] - 1] + G.KAO[i] - G.KAO[i - 1];
            for (j = G.KAO[i - 1]; j < G.KAO[i]; j++)
            {
                FO[l] = S[G.FO[j]];
                PArray[l] = G.PArray[j];
                ++l;
            }
            break;
        }
    }
    Result.VertexCount = G.VertexCount-1;
    Result.EdgeCount = l / 2;
    for (i = 0; i <= G.VertexCount; i++) {
        Result.KAO[i] = KAO[i];
    }
    for (i = 0; i < G.EdgeCount * 2; i++) {
        Result.FO[i] = FO[i];
        Result.PArray[i] = PArray[i];
    }
    Result.Targets[0] = 0;
    for (int i{ 1 }; i < v; i++) {
        Result.Targets[i] = G.Targets[i];
    }
    for (int i{ v + 1 }; i < G.VertexCount + 1; i++) {
        Result.Targets[i - 1] = G.Targets[i];
    }
    if (G.Targets[u] == 1 || G.Targets[v] == 1) Result.Targets[u] = 1;

    return Result;
}

int FindIndexEdgeForVertex(Graph G, int u, int v)
{
    for (int i{ G.KAO[u - 1] }; i < G.KAO[u]; i++) {
        if (G.FO[i] == v) {
            return i;
        }
    }
    return G.EdgeCount * 2;
}

Graph KParallelSeriesTransformation(Graph G, double& p) //В delphi процедура по сути функция void, но не в этом случае
{
    bool b = false; // обозвать бы это адекватно
    int  u = 0, v = 0, w, e1, e2;
    double p1, p2, p3; //что за нумерация??? и это по-человечески. e4 и e3 в коде вообще не использовались
    p = 1;
    /*for (int i{ 1 }; i < G.VertexCount + 1; i++) {
        if (G.KAO[i] - G.KAO[i - 1] == 1) {
            if (u == 0) u = i;
            else {
                v = i;
                if (G.Targets[i] == 1) {
                    p = p * G.PArray[G.KAO[i - 1]];
                    break;
                }
            }
            if (G.Targets[i] = 1) p = p * G.PArray[G.KAO[i - 1]];
        }
    }
    if (u != 0)
    {
        w = G.FO[G.KAO[u - 1]];
        e1 = G.KAO[u - 1];
        e2 = FindIndexEdgeForVertex(G, w, u);
        if (u < w) {
            for (int i{ u }; i < w - 1; i++)
                G.KAO[i] = G.KAO[i + 1] - 1;

            for (int i{ w - 1 }; i < G.VertexCount; i++)
                G.KAO[i] = G.KAO[i + 1] - 2;

            for (int i{ e1 }; i < e2 - 1; i++)
                G.FO[i] = G.FO[i + 1];

            for (int i{ e2 - 1 }; i < G.EdgeCount * 2 - 1; i++)
                G.FO[i] = G.FO[i + 2];

            for (int i{ 0 }; i < G.EdgeCount * 2 - 1; i++) {
                if (G.FO[i] > u) G.FO[i]--;
            }

            for (int i{ e1 }; i < e2; i++)
                G.PArray[i] = G.PArray[i + 1];

            for (int i{ e2 }; i < G.EdgeCount * 2 - 1; i++)
                G.PArray[i] = G.PArray[i + 2];
        }
        else {
            for (int i{ w }; i < u; i++)
                G.KAO[i]--;

            for (int i{ u }; i < G.VertexCount; i++)
                G.KAO[i] = G.KAO[i + 1] - 2;

            for (int i{ e2 }; i < e1; i++)
                G.FO[i] = G.FO[i + 1];

            for (int i{ e1 }; i < G.EdgeCount * 2 - 1; i++)
                G.FO[i] = G.FO[i + 2];

            for (int i{ 0 }; i < G.EdgeCount * 2 - 1; i++) {
                if (G.FO[i] > u) G.FO[i]--;
            }
            for (int i{ e2 }; i < e1 - 1; i++)
                G.PArray[i] = G.PArray[i + 1];
            for (int i{ e1 - 1 }; i < G.EdgeCount * 2 - 1; i++)
                G.PArray[i] = G.PArray[i + 2];
        }
        if (G.Targets[u] == 1) G.Targets[w] = 1;
        for (int i{ u }; i < G.VertexCount; i++)
            G.Targets[i] = G.Targets[i + 1];
        G.VertexCount--;
        G.EdgeCount--;
    }*/
    //+ дальше идет вот это KPST без единицы. Ничего не понимаю.
    //KPST?
    int i;
    for (i = 1; i < G.VertexCount + 1; i++) {
        if ((G.KAO[i] - G.KAO[i - 1] == 2) && ((G.Targets[i] == 0) ||
            ((G.Targets[i] == 1) && (G.Targets[G.FO[G.KAO[i - 1]]] == 1) &&
                (G.Targets[G.FO[G.KAO[i] - 1]] == 1)))) {
            b = true;
            v = i;
            u = G.FO[G.KAO[i - 1]];
            w = G.FO[G.KAO[i] - 1];
            p1 = G.PArray[G.KAO[i - 1]];
            p2 = G.PArray[G.KAO[i] - 1];
            break;
        }
    }
    while (b) {
        if (G.Targets[i] == 1) {
            p *= (p1 + p2 - p1 * p2); // выдает "неинициализированна" => не зашло в длинный if ЭТО ДЕЛЬТА
            p3 = p1 * p2 / (p1 + p2 - p1 * p2); //rc
        }
        else p3 = p1 * p2;
        G = Transformation(G, u, v, w, p3);
        b = false;
        for (i = 1; i < G.VertexCount + 1; i++) {
            if (G.KAO[i] - G.KAO[i - 1] == 2 && ((G.Targets[i] == 0) ||
                ((G.Targets[i] == 1) && (G.Targets[G.FO[G.KAO[i - 1]]] == 1) &&
                    (G.Targets[G.FO[G.KAO[i] - 1]] == 1)))) {
                b = true;
                v = i;
                u = G.FO[G.KAO[i - 1]];
                w = G.FO[G.KAO[i] - 1];
                p1 = G.PArray[G.KAO[i - 1]];
                p2 = G.PArray[G.KAO[i] - 1];
                break;
            }
        }
    }
    
    return G;
}

Graph Transformation(Graph G, int u, int v, int w, double &p) /*Тут был недоумевающий комментарий. Сами гадайте, что я там не поняла*/
{
    int l = 0;  
    //int targetsLength = G.VertexCount, kaoLength = G.VertexCount, foLength, pArrayLength; //размеры для Result. Сделаем пока так, чтобы не запутаться с G
    int Numbers[NV+1];
    Graph Result;
    Result.VertexCount = G.VertexCount - 1;
    for (int i{ 1 }; i < v; i++) {
        Result.Targets[i] = G.Targets[i];
    }
    for (int i{ v }; i < Result.VertexCount+1; i++) {
        Result.Targets[i] = G.Targets[i + 1];
    }
     //SetLength(Result.KAO, length(G.KAO)-1);
    Result.KAO[0] = 0;
    for (int i{ 1 }; i < v;i++) {
        Numbers[i] = i;
    }
    Numbers[v] = 0;
    for (int i{ v + 1 }; i < G.VertexCount + 1; i++) {
        Numbers[i] = i - 1;
    }
    //С этого момента начитается хоррор
    /*Загадка человечества - много дублирования
      Довольно запутанная ситуация из-за этого, стоит проверить раз 10*/
    int SE = FindIndexEdgeForVertex(G,u,w);
    bool bol = true;
    if (SE < G.EdgeCount * 2) {
        p = 1 - (1 - G.PArray[SE]) * (1 - p); //Формула? Какая?
        //SetLength(Result.FO, length(G.FO)-4); 
        Result.EdgeCount = G.EdgeCount - 2;
        for (int i{ 1 }; i < G.VertexCount + 1;i++) {
            if (i == u || i == w) {
                Result.KAO[Numbers[i]] = Result.KAO[Numbers[i] - 1];
                for (int j{ G.KAO[i - 1] }; j < G.KAO[i]; j++) { 
                    if (G.FO[j] == w || G.FO[j] == u) {
                        Result.FO[l] = Numbers[G.FO[j]];
                        Result.PArray[l] = p;
                        Result.KAO[Numbers[i]]++;
                        l++;
                    }
                    else if (G.FO[j] != v) {
                        Result.FO[l] = Numbers[G.FO[j]];
                        Result.PArray[l] = G.PArray[j];;
                        Result.KAO[Numbers[i]]++;
                        l++;
                    }
                }
            }
            else if (i != v) {
                Result.KAO[Numbers[i]] = Result.KAO[Numbers[i] - 1];
                for (int j{ G.KAO[i - 1] }; j < G.KAO[i]; j++) {
                    Result.FO[l] = Numbers[G.FO[j]];
                    Result.PArray[l] = G.PArray[j];;
                    Result.KAO[Numbers[i]]++;
                    l++;
                }
            }
        }
    }
    else {
        /*foLength = G.EdgeCount * 2 - 2;
        pArrayLength = G.EdgeCount * 2 - 2;*/
        //SetLength(Result.FO, length(G.FO)-2);
        Result.EdgeCount = G.EdgeCount - 1;
        for (int i{ 1 }; i < G.VertexCount + 1; i++) {
            if (i == u) {
                Result.KAO[Numbers[i]] = Result.KAO[Numbers[i] - 1];
                Result.FO[l] = Numbers[w];
                Result.PArray[l] = p;
                Result.KAO[Numbers[i]]++;
                l++;
                for (int j{ G.KAO[i - 1] }; j < G.KAO[i]; j++) { //bpv
                    if (G.FO[j] != v) {
                        Result.FO[l] = Numbers[G.FO[j]];
                        Result.PArray[l] = G.PArray[j];
                        Result.KAO[Numbers[i]]++;
                        l++;
                    }
                }
            }
            else if (i == w) {
                Result.KAO[Numbers[i]] = Result.KAO[Numbers[i] - 1];
                Result.FO[l] = Numbers[u];
                Result.PArray[l] = p;
                Result.KAO[Numbers[i]]++;
                l++;
                for (int j{ G.KAO[i - 1] }; j < G.KAO[i]; j++) {//bpv
                    if (G.FO[j] != v) {
                        Result.FO[l] = Numbers[G.FO[j]];
                        Result.PArray[l] = G.PArray[j];
                        Result.KAO[Numbers[i]]++;
                        l++;
                    }
                }
            }
            else if (i != v) {
                Result.KAO[Numbers[i]] = Result.KAO[Numbers[i] - 1];
                for (int j{ G.KAO[i - 1] }; j < G.KAO[i]; j++) {
                    Result.FO[l] = Numbers[G.FO[j]];
                    Result.PArray[l] = G.PArray[j];
                    Result.KAO[Numbers[i]]++;
                    l++;
                }
            }
        }
    }
    //Result.VertexCount = kaoLength-1;
    //Result.EdgeCount = foLength / 2;
    return Result;
}

bool KComponent(Graph &G)
{
    int sum, sum1 = 0, sumAll, l, LB;
    bool result;
    int Spot[NV + 1], B[NV + 1], A[NV + 1];
    for (int i{ 1 }; i < G.VertexCount+1; i++) {
        if (G.Targets[i] == 1)
            sum1++;
    }
    if (sum1 == 1) {
        G.VertexCount = 2;
        G.KAO[0] = 0;
        G.KAO[1] = 0;
        G.EdgeCount = 0;
        G.Targets[0] = 0;
        G.Targets[1] = 0;
        result = true;
    }
    else {
        sum = 1;
        sumAll = 1;
        //SetLength(A,1)
        for (int i{ 1 }; i < G.VertexCount + 1; i++) {
            if (G.Targets[i] == 1) {
                A[0] = i;
                break;
            }
        }
        l = 1;
        //SetLength(Spot, length(G.KAO))
        for (int i{ 1 }; i < G.VertexCount + 1; i++) {
            Spot[i] = 2;
        }
        Spot[A[0]] = 1;
        while (l > 0) {
            //SetLength(B,0);
            LB = 0;
            for (int i{ 0 };i < l;i++) {
                for (int j{ G.KAO[A[i] - 1] }; j < G.KAO[A[i]]; j++) {
                    if (Spot[G.FO[j]] == 2) {
                        LB++;
                        B[LB - 1] = G.FO[j];
                        Spot[G.FO[j]] = 1;
                        if (G.Targets[G.FO[j]] == 1) sum++;
                        sumAll++;
                    }
                }
            }
            l = LB;
            if (l > 0) {
                for (int i{ 0 }; i < l;i++) {
                    A[i] = B[i];
                }
            }
        }

    }
    if (sum == sum1) result = true;
    else result = false;
    //if (result = true && sumAll < G.VertexCount + 1) G = InducedKGraph(G, Spot, 1);

    return result;
    //if (Result=true) and (SumAll<(length(G.KAO)-1)) then G:=InducedKGraph(G,Spot,1);?????
}

Graph InducedKGraph(Graph G, int list[], int number)
{
    Graph Result;
    int S[NV + 1], FO[NE * 2], KAO[NV + 1], l;
    double PArray[NE * 2];
    auto boolka{ [](Graph G, int list[], int number, int i, int j) {
            if (number == 1) return true;
            else {
                if (list[i] == 0 && list[G.FO[j]] == 0) return false;
                else return true;
            }
        } };
    int j = 0;
    for (int i{ 1 }; i < G.VertexCount + 1; i++) {
        if (list[i] == number || list[i] == 0) {
            j++;
            S[i] = j;
        }
        else S[i] = 0;
    }
    //SetLength(FO,length(G.FO));
    //SetLength(Parray,length(G.FO));
    KAO[0] = 0;
    Result.VertexCount = j + 1;
    l = 0;
    for (int i{ 1 }; i < G.VertexCount + 1; i++) {
        if (S[i] != 0) {
            KAO[S[i]] = KAO[S[i] - 1];
            for (int j{ G.KAO[i - 1] }; i < G.KAO[i] - 1; i++) {
                if (S[G.FO[j]] != 0 && boolka(G, list, number, i, j)) {
                    KAO[S[i]]++;
                    FO[l] = S[G.FO[j]];
                    PArray[l] = G.PArray[j];
                    l++;
                }
            }
        }
    }
    //SetLength(FO,l);
    //SetLength(Parray, l);
    Result.EdgeCount = l / 2;
    for (int i{ 1 }; i < Result.VertexCount + 1; i++) {
        Result.KAO[i] = KAO[i];
    }
    for (int i{ 0 }; i < G.EdgeCount * 2; i++) {
        Result.FO[i] = FO[i];
        Result.PArray[i] = PArray[i];
    }
    Result.Targets[0] = 0;
    for (int i{ 1 }; i < Result.VertexCount + 1; i++) {
        if (S[i] > 0) Result.Targets[S[i]] = G.Targets[i];
    }
    return Result;
}



