#include <algorithm>
#include <float.h>
#include <assert.h>
#include <stdio.h>
#include "rect_core.hpp"


static bool isSimilarRects(const Rect_T& r1, 
                           const Rect_T& r2,
                           double eps);

static int clusterRects(const std::vector<Rect_T>& rects,
                        std::vector<int>& labels, 
                        double eps);

Rect_T intersect(const Rect_T &r1, const Rect_T &r2)
{
    Rect_T res = {0,0,0,0};
    res.x = MAX_T(r1.x, r2.x);
    res.y = MAX_T(r1.y, r2.y);
    res.w = MIN_T(r1.x+r1.w, r2.x+r2.w)-res.x;
    res.h = MIN_T(r1.y+r1.h, r2.y+r2.h)-res.y;    
    return res;
}

void groupRects(std::vector<Rect_T>& rectList, 
                int groupThreshold, double eps, 
                std::vector<int>* weights, 
                std::vector<float>* pConfs)
{
    if (groupThreshold<=0 || rectList.empty())
    {
        if (weights)
        {
            int sz = (int)rectList.size();
            weights->resize(sz);
            for(int i=0; i<sz; i++ )
                (*weights)[i] = 1;
        }
        return;
    }
    
    if (0!=pConfs && pConfs->empty())
        pConfs = 0;

    std::vector<int> labels;
    int nclasses = clusterRects(rectList, labels, eps);

    std::vector<Rect_T> rrects(nclasses);
    std::vector<int> rweights(nclasses, 0);
    std::vector<float> rconfs(nclasses, 0);
    int i, j, nlabels = (int)labels.size();
    for( i = 0; i < nlabels; i++ )
    {
        int cls = labels[i];
        rrects[cls].x += rectList[i].x;
        rrects[cls].y += rectList[i].y;
        rrects[cls].w += rectList[i].w;
        rrects[cls].h += rectList[i].h;
        rweights[cls]++;
        if (0 != pConfs)
            rconfs[cls] += (*pConfs)[i]; 
    }
    
    for( i = 0; i < nclasses; i++ )
    {
        Rect_T r = rrects[i];
        float s = 1.f/rweights[i];
        rrects[i].x = static_cast<int>(r.x*s);
        rrects[i].y = static_cast<int>(r.y*s);
        rrects[i].w = static_cast<int>(r.w*s);
        rrects[i].h = static_cast<int>(r.h*s);
        if (0 != pConfs)
            rconfs[i] *= s; 
    }

    rectList.clear();
    if (0 != weights)
        weights->clear();
    if (0 != pConfs)
        pConfs->clear();
    
    for( i = 0; i < nclasses; i++ )
    {
        Rect_T r1 = rrects[i];
        int n1 = rweights[i];
        if (n1 <= groupThreshold )
            continue;
        
        //filter out small rects inside large rects
        for( j = 0; j < nclasses; j++ )
        {
            int n2 = rweights[j];            
            if( j == i || n2 <= groupThreshold )
                continue;

            Rect_T r2 = rrects[j];            
            int dx = static_cast<int>(r2.w * eps);
            int dy = static_cast<int>(r2.h * eps);

            if (r1.x >= r2.x - dx &&
                r1.y >= r2.y - dy &&
                r1.x + r1.w <= r2.x + r2.w + dx &&
                r1.y + r1.h <= r2.y + r2.h + dy &&
                (n2 > MAX_T(3, n1) || n1 < 3) )
                break;
        }
        
        if( j == nclasses )
        {
            rectList.push_back(r1);
            if (0 != weights)
                weights->push_back(n1);
            if (0 != pConfs)
                pConfs->push_back(rconfs[i]);
        }
    }
}

static bool isSimilarRects(const Rect_T& r1, 
                           const Rect_T& r2,
                           double eps)
{
    double delta = eps*(MIN_T(r1.w, r2.w) + 
                        MIN_T(r1.h, r2.h))*0.5;
#if defined(__APPLE__) 
    if (abs(r1.x-r2.x) <= delta &&
        abs(r1.y-r2.y) <= delta &&
        abs(r1.x+r1.w-r2.x-r2.w) <= delta &&
        abs(r1.y+r1.h-r2.y-r2.h) <= delta)
        return true;
#else
    if (std::abs(r1.x-r2.x) <= delta &&
        std::abs(r1.y-r2.y) <= delta &&
        std::abs(r1.x+r1.w-r2.x-r2.w) <= delta &&
        std::abs(r1.y+r1.h-r2.y-r2.h) <= delta)
        return true;
#endif
    return false;
}


static int clusterRects(const std::vector<Rect_T>& rects,
                        std::vector<int>& labels, double eps)
{
    int i, j, N = (int)rects.size();
    const Rect_T* pRect = &rects[0];

    const int PARENT=0;
    const int RANK=1;

    std::vector<int> _nodes(N*2);
    int (*nodes)[2] = (int(*)[2])&_nodes[0];

    // The first O(N) pass: create N single-vertex trees
    for (i=0; i<N; i++)
    {
        nodes[i][PARENT]=-1;
        nodes[i][RANK] = 0;
    }

    // The main O(N^2) pass: merge connected components
    for (i=0; i<N; i++ )
    {
        int root = i;

        // find root
        while( nodes[root][PARENT] >= 0 )
            root = nodes[root][PARENT];

        for( j = 0; j < N; j++ )
        {
            if (i==j || !isSimilarRects(pRect[i], 
                                        pRect[j], eps))
                continue;
            int root2 = j;

            while( nodes[root2][PARENT] >= 0 )
                root2 = nodes[root2][PARENT];

            if( root2 != root )
            {
                // unite both trees
                int rank = nodes[root][RANK];
                int rank2 = nodes[root2][RANK];
                if( rank > rank2 )
                    nodes[root2][PARENT] = root;
                else
                {
                    nodes[root][PARENT] = root2;
                    nodes[root2][RANK] += rank == rank2;
                    root = root2;
                }
                assert( nodes[root][PARENT] < 0 );

                int k=j, parent;

                // compress the path from node2 to root
                while( (parent = nodes[k][PARENT]) >= 0 )
                {
                    nodes[k][PARENT] = root;
                    k = parent;
                }

                // compress the path from node to root
                k = i;
                while( (parent = nodes[k][PARENT]) >= 0 )
                {
                    nodes[k][PARENT] = root;
                    k = parent;
                }
            }
        }
    }

    // Final O(N) pass: enumerate classes
    labels.resize(N);
    int nclasses = 0;
    
    for( i = 0; i < N; i++ )
    {
        int root = i;
        while( nodes[root][PARENT] >= 0 )
            root = nodes[root][PARENT];
        // re-use the rank as the class label
        if( nodes[root][RANK] >= 0 )
            nodes[root][RANK] = ~nclasses++;
        labels[i] = ~nodes[root][RANK];
    }
    
    return nclasses;
}
