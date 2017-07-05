#ifndef __FEATURE_HPP__
#define __FEATURE_HPP__

#include <assert.h>
#include "integral.hpp"

typedef struct _tagFeaType
{
    int typeID;
    int x;
    int y;
    int ew;                    
    int eh;                    
    int ch;                     
}FeaType;

template<class T>
T Fea_Haar(const FeaType &F,
           const T *II,
           const T *II_T,
           const int W, const int H)
{
#define AREA(l1,l2,w) ((l1)[0]+(l2)[w]-(l1)[w]-(l2)[0])
    
    const T *line1;
    const T *line2;
    T white = 0, black = 0;
    switch (F.typeID) {
    case 1:
        {
            // ---+++
            // ---+++           
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + F.eh*IntegralImage_Pitch(W);
            black = AREA(line1, line2, F.ew);

            line1 += F.ew;
            line2 += F.ew;
            white = AREA(line1, line2, F.ew);            
            return (white - black);
        }
    case 2: 
        {
            // ---
            // ---
            // +++
            // +++
            int pitch = F.eh*IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + pitch;
            black = AREA(line1, line2, F.ew);

            line1 = line2 + pitch;
            white = AREA(line2, line1, F.ew);
            return (white - black);
        }
    case 3:
        {
            // ---+++---
            // ---+++---
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + F.eh*IntegralImage_Pitch(W);
            black = AREA(line1, line2, F.ew*3);

            line1 += F.ew;
            line2 += F.ew;
            white = AREA(line1, line2, F.ew);
            return (3*white - black);
        }
    case 4:
        {
            // ---
            // ---
            // +++
            // +++
            // ---
            // ---
            int pitch = F.eh*IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + 3*pitch;
            black = AREA(line1, line2, F.ew);

            line1 += pitch;
            line2 -= pitch;
            white = AREA(line1, line2, F.ew);
            return (3*white - black);
        }
    case 5:
        {
            // +++---
            // +++---
            // ---+++
            // ---+++
            int pitch;
            const T* line3;

            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            pitch = F.eh * IntegralImage_Pitch(W);
            line2 = line1 + pitch;
            white = AREA(line1, line2, F.ew);

            line3 = line2 + pitch;
            black = AREA(line1, line3, (F.ew<<1));

            line2 += F.ew;
            line3 += F.ew;
            white += AREA(line2, line3, F.ew);
            return (2*white - black);
        }
    case 6:
        {
            // ---++++++---
            // ---++++++---
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + F.eh*IntegralImage_Pitch(W);
            black = AREA(line1, line2, (F.ew<<2));

            line1 += F.ew;
            line2 += F.ew;
            white = AREA(line1, line2, (F.ew<<1));
            return (2*white - black);
        }
    case 7:
        {
            // ---
            // ---
            // +++
            // +++
            // +++
            // +++
            // ---
            // ---
            int pitch = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + (pitch<<2);
            black = AREA(line1, line2, F.ew);

            line1 += pitch;
            line2 -= pitch;
            white = AREA(line1, line2, F.ew);
            return (2 * white - black);
        }
    case 8:
        {
            // ---------
            // ---------
            // ---+++---
            // ---+++---
            // ---------
            // ---------
            int pitch = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(II, W, F.y) + F.x;
            line2 = line1 + 3*pitch;
            black = AREA(line1, line2, F.ew*3);

            line1 += pitch + F.ew;
            line2 -= pitch - F.ew;
            white = AREA(line1, line2, F.ew);
            return (9 * white - black);
        }
    }
    return 0;

#undef AREA
}

#ifdef ENABLE_NEON
int Fea_MBLBP_4P_Neon(const FeaType &F,
                      const unsigned int *II,
                      const int W, const int H,
                      unsigned int *respond);
#endif
    
template<class T>
T Fea_MBLBP(const FeaType &F,
            const T *II,
            const int W, const int H)
{    
    const T *pLine;
    int pitch = F.eh*IntegralImage_Pitch(W);
    
    T respond = 0;
    T p01, p02, p03, p04, p05, p06, p07, p08;
    T p09, p10, p11, p12, p13, p14, p15, p16;
    T val0, val1, val2, val3, cent, val4, val5, val6, val7;

    //assert(F.x >= 0);
    //assert(F.y >= 0);
    //assert(F.x + 3*F.ew <= W);
    //assert(F.y + 3*F.eh <= H);

    pLine = IntegralImage_Line(II, W, F.y) + F.x;
    int w1=F.ew;
    int w2=w1<<1;
    int w3=w2+w1;
    p01 = pLine[0];
    p02 = pLine[w1];
    p03 = pLine[w2];
    p04 = pLine[w3];

    pLine +=  pitch;
    p05 = pLine[0];
    p06 = pLine[w1];
    p07 = pLine[w2];
    p08 = pLine[w3];

    pLine +=  pitch;
    p09 = pLine[0];
    p10 = pLine[w1];
    p11 = pLine[w2];
    p12 = pLine[w3];

    pLine +=  pitch;
    p13 = pLine[0];
    p14 = pLine[w1];
    p15 = pLine[w2];
    p16 = pLine[w3];

    val0 = p06+p01-p02-p05;
    val1 = p07+p02-p03-p06;
    val2 = p08+p03-p04-p07;
    val3 = p10+p05-p06-p09;
    cent = p11+p06-p07-p10;
    val4 = p12+p07-p08-p11;
    val5 = p14+p09-p10-p13;
    val6 = p15+p10-p11-p14;
    val7 = p16+p11-p12-p15;
    
    if (val0>cent) 
        respond = respond+1;
    if (val1>cent) 
        respond = respond+2;
    if (val2>cent)
        respond = respond+4;
    if (val3>cent) 
        respond = respond+8;
    if (val4>cent)
        respond = respond+16;
    if (val5>cent)
        respond = respond+32;
    if (val6>cent)
        respond = respond+64;
    if (val7>cent) 
        respond = respond+128;
    
    return respond;
}

template<class T>
T Fea_Haar_backup(const FeaType &F,
           const T *II,
           const T *II_T,
           const int W, const int H)
{
#define AREA_TILTED(y,x) IntegralImage_Value(II_T, W, x, y)
#define AREA(y,x) IntegralImage_Value(II, W, x, y)
#define AREA_fast(l1,l2,w) ((l1)[0]+(l2)[w]-(l1)[w]-(l2)[0])
    
    const T * const ii = IntegralImage_TopLeft(II, W);
    const T *line1;
    const T *line2;
    T white = 0, black = 0;
    const int x = F.x - 1;
    const int y = F.y - 1;
    assert(x >= -1);
    assert(y >= -1);

    switch (F.typeID) {
    case 1:
        {
            // ---+++
            // ---+++           
            assert(x + 2*F.ew < W);
            assert(y + F.eh < H);
            
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + F.eh * IntegralImage_Pitch(W);
            black = AREA_fast(line1, line2, F.ew);
            line1 += F.ew;
            line2 += F.ew;
            white = AREA_fast(line1, line2, F.ew);
            
            return (white - black);
        }
    case 2: 
        {
            // ---
            // ---
            // +++
            // +++
            assert(x + F.ew < W);
            assert(y + 2*F.eh < H);
            
            int s = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + s;
            black = AREA_fast(line1, line2, F.ew);
            line1 = line2 + s;
            white = AREA_fast(line2, line1, F.ew);
            return (white - black);
        }
    case 3:
        {
            // ---+++---
            // ---+++---
            assert(x + 3*F.ew < W);
            assert(y + F.eh < H);
            
            int w = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + w;
            w = 3 * F.ew;
            black = AREA_fast(line1, line2, w);
            line1 += F.ew;
            line2 += F.ew;
            white = AREA_fast(line1, line2, F.ew);
            return (3*white - black);
        }
    case 4:
        {
            // ---
            // ---
            // +++
            // +++
            // ---
            // ---
            assert(x + F.ew < W);
            assert(y + 3*F.eh < H);
            
            int s = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + 3 * s;
            black = AREA_fast(line1, line2, F.ew);
            line1 += s;
            line2 -= s;
            white = AREA_fast(line1, line2, F.ew);
            return (3*white - black);
        }
    case 5:
        {
            // +++---
            // +++---
            // ---+++
            // ---+++
            assert(x + 2*F.ew < W);
            assert(y + 2*F.eh < H);
            int w, s;
            const T* line3;
            line1 = IntegralImage_Line(ii, W, y) + x;
            w = 2 * F.ew;
            s = F.eh * IntegralImage_Pitch(W);
            line2 = line1 + s;
            line3 = line1 + (s<<1);
            black = AREA_fast(line1, line3, w);
            white = AREA_fast(line1, line2, F.ew);
            line2 += F.ew;
            line3 += F.ew;
            white += AREA_fast(line2, line3, F.ew);
            return (2*white - black);
        }
    case 6:
        {
            // ---++++++---
            // ---++++++---
            assert(x+4*F.ew < W);
            assert(y+F.eh < H);
            
            int w = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + w;
            w = 4 * F.ew;
            black = AREA_fast(line1, line2, w);
            line1 += F.ew;
            line2 += F.ew;
            w >>= 1;
            white = AREA_fast(line1, line2, w);
            return (2*white - black);
        }
    case 7:
        {
            // ---
            // ---
            // +++
            // +++
            // +++
            // +++
            // ---
            // ---
            assert(x+F.ew < W);
            assert(y+4*F.eh < H);
            
            int s = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + (s<<2);
            black = AREA_fast(line1, line2, F.ew);
            line1 += s;
            line2 -= s;
            white = AREA_fast(line1, line2, F.ew);
            return (2 * white - black);
        }
    case 8:
        {
            // ---------
            // ---------
            // ---+++---
            // ---+++---
            // ---------
            // ---------
            assert(x + 3*F.ew < W);
            assert(y + 3*F.eh < H);
            
            int w = 3 * F.ew;
            int s = F.eh * IntegralImage_Pitch(W);
            line1 = IntegralImage_Line(ii, W, y) + x;
            line2 = line1 + 3 * s;
            black = AREA_fast(line1, line2, w);
            line1 += s + F.ew;
            line2 = line2 - s + F.ew;
            white = AREA_fast(line1, line2, F.ew);
            return (9 * white - black);
        }
        /////////////////////////////////////////
        //following are tilted haar-like features
    case 9:      
        {
            //+++---
            //+++---         
            assert((x-F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew*2+1<W)&&(y+F.ew*2+F.eh<H));
            
            white = AREA_TILTED(y,x+1)  
                + AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1) 
                - AREA_TILTED(y+F.ew,x+F.ew+1)  
                - AREA_TILTED(y+F.eh,x-F.eh+1); 
            
            black = AREA_TILTED(y+F.ew,x+F.ew+1)
                + AREA_TILTED(y+F.eh+F.ew*2,x-F.eh+F.ew*2+1) 
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1) 
                - AREA_TILTED(y+F.ew*2,x+F.ew*2+1);
            return (white-black);
        }
    case 10:
        {
            //+++
            //+++
            //---
            //---
            assert((x-2*F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew+1<W)&&(y+F.ew+F.eh*2<H));

            white = AREA_TILTED(y,x+1) 
                + AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1)
                - AREA_TILTED(y+F.ew,x+F.ew+1) 
                - AREA_TILTED(y+F.eh,x-F.eh+1); 
            
            black = AREA_TILTED(y+F.eh,x-F.eh+1)  
                +  AREA_TILTED(y+F.eh*2+F.ew,x-F.eh*2+F.ew+1)
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1)
                - AREA_TILTED(y+F.eh*2,x-F.eh*2+1);; 
            return (white-black);
        }
    case 11:
        {
            //+++---+++
            //+++---+++
            assert((x-F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew*3+1<W)&&(y+F.ew*3+F.eh<H));
        
            white = AREA_TILTED(y,x+1)  
                + AREA_TILTED(y+F.eh+F.ew*3,x-F.eh+F.ew*3+1)
                - AREA_TILTED(y+F.ew*3,x+F.ew*3+1)   
                - AREA_TILTED(y+F.eh,x-F.eh+1); 
            
            black = AREA_TILTED(y+F.ew,x+F.ew+1)
                + AREA_TILTED(y+F.eh+F.ew*2,x-F.eh+F.ew*2+1) 
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1) 
                - AREA_TILTED(y+F.ew*2,x+F.ew*2+1);
            return (white-3*black);
        }
    case 12:
        {
            //+++-----+++
            //+++-----+++
            assert((x-F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew*4+1<W)&&(y+F.ew*4+F.eh<H));
        
            white = AREA_TILTED(y,x+1)  
                + AREA_TILTED(y+F.eh+F.ew*4,x-F.eh+F.ew*4+1) 
                - AREA_TILTED(y+F.ew*4,x+F.ew*4+1)  
                - AREA_TILTED(y+F.eh,x-F.eh+1); 
            
            black = AREA_TILTED(y+F.ew,x+F.ew+1)
                + AREA_TILTED(y+F.eh+F.ew*3,x-F.eh+F.ew*3+1)
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1) 
                - AREA_TILTED(y+F.ew*3,x+F.ew*3+1);
            return (white-2*black);
        }
    case 13:
        {
            //+++
            //+++
            //---
            //---
            //+++
            //+++
            assert((x-3*F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew+1<W)&&(y+F.ew+F.eh*3<H));

            white = AREA_TILTED(y,x+1) 
                + AREA_TILTED(y+F.eh*3+F.ew,x-F.eh*3+F.ew+1)
                - AREA_TILTED(y+F.ew,x+F.ew+1) 
                - AREA_TILTED(y+F.eh*3,x-F.eh*3+1); 
            
            black = AREA_TILTED(y+F.eh,x-F.eh+1)  
                +  AREA_TILTED(y+F.eh*2+F.ew,x-F.eh*2+F.ew+1)
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1)
                - AREA_TILTED(y+F.eh*2,x-F.eh*2+1);; 
            return (white-3*black);
        }
    case 14:
        {
            //+++
            //+++
            //---
            //---
            //---
            //---
            //+++
            //+++
            assert((x-4*F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew+1<W)&&(y+F.ew+F.eh*4<H));

            white = AREA_TILTED(y,x+1) 
                + AREA_TILTED(y+F.eh*4+F.ew,x-F.eh*4+F.ew+1)
                - AREA_TILTED(y+F.ew,x+F.ew+1) 
                - AREA_TILTED(y+F.eh*4,x-F.eh*4+1);
            
            black = AREA_TILTED(y+F.eh,x-F.eh+1)  
                +  AREA_TILTED(y+F.eh*3+F.ew,x-F.eh*3+F.ew+1)
                - AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1)
                - AREA_TILTED(y+F.eh*3,x-F.eh*3+1);; 
            return (white-3*black);
        }
    case 15:
        {
            //+++++++++
            //+++++++++
            //+++---+++
            //+++---+++
            //+++++++++
            //+++++++++
            assert((x-3*F.eh+1>=-1)&&(y>=-1));
            assert((x+F.ew*3+1<W)&&(y+F.ew*3+F.eh*3<H));

            white = AREA_TILTED(y,x+1) 
                + AREA_TILTED(y+F.eh*3+F.ew*3,x-F.eh*3+F.ew*3+1) 
                - AREA_TILTED(y+F.ew*3,x+F.ew*3+1) 
                - AREA_TILTED(y+F.eh*3,x-F.eh*3+1);
        
            black = AREA_TILTED(y+F.eh+F.ew,x-F.eh+F.ew+1) 
                + AREA_TILTED(y+F.eh*2+F.ew*2,x-F.eh*2+F.ew*2+1) 
                - AREA_TILTED(y+F.eh+F.ew*2,x-F.eh+F.ew*2+1)
                - AREA_TILTED(y+F.eh*2+F.ew,x-F.eh*2+F.ew+1);
            return (white-9*black);
        }        
    }

    return 0;

#undef AREA
#undef AREA_fast
}

#endif
