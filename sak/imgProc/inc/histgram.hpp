#ifndef __HISTGRAM_HPP__
#define __HISTGRAM_HPP__
#include <string.h>
#include <math.h>

template<typename T1, typename T2> 
int calcHist2D(T1 *src0, T1 *src1, 
               int width, int height, int step, 
               T2 *pHist, int *bins, float **ranges);

template<typename T1, typename T2> 
int calcHist2D_epan(T1 *src0, T1 *src1, 
                    int width, int height, int step, 
                    T2 *pHist, int *bins, float **ranges);

template<typename T1, typename T2> 
int calcHist3D(T1 *src0, T1 *src1, T1 *src2,
               int width, int height, int step, 
               T2 *pHist, int *bins, float **ranges);

template<typename T>
int normHist_MinMax(T *pHist, int len, T minV, T maxV);

template<typename T>
int normHist_Norm(T *pHist, int len);

template<typename T1, typename T2>
int calcBackProj2D(T1 *src0, T1 *src1, 
                   int width, int height, int step, 
                   T2  *roiHist, int *bins, float **ranges,
                   unsigned char *probImg);

template<typename T1, typename T2>
int calcBackProj3D(T1 *src0, T1 *src1, T1 *src2, 
                   int width, int height, int step, 
                   T2  *roiHist, int *bins, float **ranges,
                   unsigned char *probImg);

template<typename T1, typename T2> 
int calcHist3D(T1 *src0, T1 *src1, T1 *src2,
               int width, int height, int step, 
               T2 *pHist, int *bins, float **ranges)
{
    /*
      ranges shoud be (0,255) when in unsigned char
     */

    if ((0==src0) || (0==src1) || (0==src2) || (0==pHist) 
        || (0==bins) ||(0==ranges))
        return -1;

    memset(pHist, 0, bins[0]*bins[1]*bins[2]*sizeof(*pHist));

    float cs[3] = {bins[0]/(ranges[0][1]-ranges[0][0]+1), 
                   bins[1]/(ranges[1][1]-ranges[1][0]+1),
                   bins[2]/(ranges[2][1]-ranges[2][0]+1)};
    
    for (int h=0; h<height; h++)
    {
        for (int w=0; w<width; w++)
        {
            int idx_0 = (int)(src0[w]*cs[0]);
            int idx_1 = (int)(src1[w]*cs[1]);
            int idx_2 = (int)(src2[w]*cs[2]);
            int idx = idx_0*bins[1]*bins[2]+idx_1*bins[2]+idx_2;
            pHist[idx] += 1;         
        }
        src0 += step;
        src1 += step;
        src2 += step;
    }
    return 0;
}

template<typename T1, typename T2> 
int calcHist2D_epan(T1 *src0, T1 *src1, 
                    int width, int height, int step, 
                    T2 *pHist, int *bins, float **ranges)
{
    /*
      ranges shoud be (0,255) when in unsigned char
    */
    
    if ((0==src0) || (0==src1) || (0==pHist) 
        || (0==bins) ||(0==ranges))
        return -1;
    
    memset(pHist, 0, bins[0]*bins[1]*sizeof(*pHist));
    
    float cs[2] = {bins[0]/(ranges[0][1]-ranges[0][0]+1), 
                   bins[1]/(ranges[1][1]-ranges[1][0]+1)};

    float centre_x = (width-1)/2.f;
    float centre_y = (height-1)/2.f;
    
    for (int h=0; h<height; h++)
    {
        float norm_y = (h-centre_y)/centre_y;
        norm_y = norm_y*norm_y;
        for (int w=0; w<width; w++)
        {
            float norm_x = (w-centre_x)/centre_x;
            norm_x = norm_x*norm_x;
            float d = norm_x+norm_y;
            if (d<1)
            {            
                int idx_0 = (int)(src0[w]*cs[0]);
                int idx_1 = (int)(src1[w]*cs[1]);
                int idx = idx_0*bins[1]+idx_1;
                // epanechnikov with cd=2, d=1
                pHist[idx] += (1-d)*10.f;   
            }      
        }
        src0 += step;
        src1 += step;
    }
    return 0;
}

template<typename T1, typename T2> 
int calcHist2D(T1 *src0, T1 *src1, 
               int width, int height, int step, 
               T2 *pHist, int *bins, float **ranges)
{
    /*
      ranges shoud be (0,255) when in unsigned char
     */

    if ((0==src0) || (0==src1) || (0==pHist) 
        || (0==bins) ||(0==ranges))
        return -1;

    memset(pHist, 0, bins[0]*bins[1]*sizeof(*pHist));

    float cs[2] = {bins[0]/(ranges[0][1]-ranges[0][0]+1), 
                   bins[1]/(ranges[1][1]-ranges[1][0]+1)};
    
    for (int h=0; h<height; h++)
    {
        for (int w=0; w<width; w++)
        {
            int idx_0 = (int)(src0[w]*cs[0]);
            int idx_1 = (int)(src1[w]*cs[1]);
            int idx = idx_0*bins[1]+idx_1;
            pHist[idx] += 1;         
        }
        src0 += step;
        src1 += step;
    }
    return 0;
}

template<typename T>
int normHist_MinMax(T *pHist, int len, T minV, T maxV)
{
    if (0==pHist)
        return -1;
    
    float range_c=0, range=maxV-minV;    
    T minV_c=pHist[0], maxV_c=pHist[0];
    T *ptr=0, val;    

    ptr = pHist;
    for (int l=0; l<len; l++)
    {
        val = ptr[l];        
        if (val < minV_c)
            minV_c = val;
        if (val > maxV_c)
            maxV_c = val;
    }
    
    range_c = maxV_c - minV_c;    
    ptr = pHist;
    for (int l=0; l<len; l++)
        ptr[l] = (ptr[l]-minV_c)/range_c*range + minV;  
    return 0;
}   

template<typename T>
int normHist_Norm(T *pHist, int len)
{
    if (0==pHist)
        return -1;
    float sum = 0;

    for (int l=0; l<len; l++)
        sum += pHist[l];

    for (int l=0; l<len; l++)
        pHist[l] = pHist[l]/sum;
    return 0;
}   

template<typename T1, typename T2>
int calcBackProj2D(T1 *src0, T1 *src1, 
                   int width, int height, int step, 
                   T2 *roiHist, int *bins, float **ranges,
                   unsigned char *probImg)
{
    /*
      ranges shoud be (0,255) when in unsigned char
     */
    if ((0==src0) || (0==src1) || (0==roiHist) 
        || (0==bins) || (0==probImg) || (0==ranges))
        return -1;

    float cs[2] = {bins[0]/(ranges[0][1]-ranges[0][0]+1), 
                   bins[1]/(ranges[1][1]-ranges[1][0]+1)};

    for (int h=0; h<height; h++)
    {
        for (int w=0; w<width; w++)
        {
            int idx_0 = (int)(src0[w]*cs[0]);
            int idx_1 = (int)(src1[w]*cs[1]);
            int idx = idx_0*bins[1]+idx_1;
            *(probImg++) = (int)(roiHist[idx]+0.5f);  
        }
        src0 += step;
        src1 += step;
    }
    return 0;
}

template<typename T1, typename T2>
int calcBackProj3D(T1 *src0, T1 *src1, T1 *src2,
                   int width, int height, int step, 
                   T2 *roiHist, int *bins, float **ranges,
                   unsigned char *probImg)
{
    /*
      ranges shoud be (0,255) when in unsigned char
     */
    if ((0==src0) || (0==src1) || (0==src2) || (0==roiHist) 
        || (0==bins) || (0==probImg) || (0==ranges))
        return -1;

    float cs[3] = {bins[0]/(ranges[0][1]-ranges[0][0]+1), 
                   bins[1]/(ranges[1][1]-ranges[1][0]+1),
                   bins[2]/(ranges[2][1]-ranges[2][0]+1)};

    for (int h=0; h<height; h++)
    {
        for (int w=0; w<width; w++)
        {
            int idx_0 = (int)(src0[w]*cs[0]);
            int idx_1 = (int)(src1[w]*cs[1]);
            int idx_2 = (int)(src2[w]*cs[2]);
            int idx = idx_0*bins[1]*bins[2]+idx_1*bins[2]+idx_2;
            *(probImg++) = (int)(roiHist[idx]+0.5f);  
        }
        src0 += step;
        src1 += step;
        src2 += step;
    }
    return 0;
}

#endif
