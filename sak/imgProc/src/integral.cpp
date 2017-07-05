#include "integral.hpp"

#ifdef ENABLE_NEON
#include <arm_neon.h>
  #if 0 // cumulative sum is not fast on neon

int integral_neon(unsigned char *img,
                  int step, int width, int height,
                  unsigned int *II)
{
    if ((0==img)||(0==II))
        return -1;
    int II_step = width+1;
    
    //first row
    memset(II, 0, sizeof(II[0]) * II_step);
    
    unsigned int *pDst_pre = II+1;
    unsigned int *pDst = II+II_step+1;
    uint16x8_t zeroVec = vdupq_n_u16(0);
    uint8x8_t zeroVec8 = vdup_n_u8(0);
    //Actually The 2rd row not need to add pre_row, because the 1st row is 0
    for (int h = 0; h<height; ++h)
    {
        int w=0;
        unsigned int preSum = 0;
        pDst[-1] = 0;
        // vector to carry over last prefix sum value to next chunk
        uint32x4_t preSum32x4 = vdupq_n_u32(0);
        for (; w+16<width; w+=16)
        {
            //load 16 elements each time
            uint8x16_t eles = vld1q_u8(img+w);
            uint8x8_t lowEles8 = vget_low_u8(eles);
            uint8x8_t highEles8 = vget_high_u8(eles);
            
            //do add with 3 pass for lowEles
            uint16x8_t lowEles16 = vaddl_u8(lowEles8, vext_u8(zeroVec8, lowEles8, 7));
            lowEles16 = vaddq_u16(lowEles16, vextq_u16(zeroVec, lowEles16, 6));
            lowEles16 = vaddq_u16(lowEles16, vextq_u16(zeroVec, lowEles16, 4));
            
            //do add with 3 pass for highEles
            uint16x8_t highEles16 = vaddl_u8(highEles8, vext_u8(zeroVec8, highEles8, 7));
            highEles16 = vaddq_u16(highEles16, vextq_u16(zeroVec, highEles16, 6));
            highEles16 = vaddq_u16(highEles16, vextq_u16(zeroVec, highEles16, 4));
            
            // 1st 4 elements
            uint32x4_t pre_row = vld1q_u32(pDst_pre);
            uint32x4_t res32x4 = vaddq_u32(vaddq_u32(vmovl_u16(vget_low_u16(lowEles16)), preSum32x4), pre_row);
            vst1q_u32(pDst, res32x4);
            
            // 2rd 4 elements
            pDst_pre += 4;
            pDst += 4;
            pre_row = vld1q_u32(pDst_pre);
            res32x4 = vaddq_u32(vmovl_u16(vget_high_u16(lowEles16)), preSum32x4);
            // update the pre last sum befor apply pre_row
            preSum = vgetq_lane_u32(res32x4, 3);
            preSum32x4 = vdupq_n_u32(preSum);
            res32x4 = vaddq_u32(res32x4, pre_row);
            vst1q_u32(pDst, res32x4);
            
            // 3rd 4 elements
            pDst_pre += 4;
            pDst += 4;
            pre_row = vld1q_u32(pDst_pre);
            res32x4 = vaddq_u32(vaddq_u32(vmovl_u16(vget_low_u16(highEles16)), preSum32x4),pre_row);
            vst1q_u32(pDst, res32x4);
            
            // last 4 elements
            pDst_pre += 4;
            pDst += 4;
            pre_row = vld1q_u32(pDst_pre);
            res32x4 = vaddq_u32(vmovl_u16(vget_high_u16(highEles16)), preSum32x4);
            //update the pre last sum
            preSum = vgetq_lane_u32(res32x4, 3);
            preSum32x4 = vdupq_n_u32(preSum);
            res32x4 = vaddq_u32(res32x4, pre_row);
            vst1q_u32(pDst, res32x4);
            
            pDst_pre += 4;
            pDst += 4;
        }
        // do the rest
        for ( ; w<width; ++w) {
            preSum += img[w];
            *(pDst++) = preSum + *(pDst_pre++);
        }
        pDst_pre += 1;
        pDst += 1;
        img += step;
    }
    return 0;
}

  #else // only sum pre line with neon

int integral_neon(unsigned char *img,
             int step, int width, int height,
             unsigned int *II)
{
    if ((0==img)||(0==II))
        return -1;
    
    int II_step = width+1;
    memset(II, 0, (II_step)*sizeof(unsigned int));
    
    unsigned char *pSrc = img;
    unsigned int *pDst_pre=II+1;
    unsigned int *pDst = pDst_pre + II_step;
    unsigned int sum = 0;
    //first line
    {
        pDst[-1] = 0;
        for (int w=0; w<width; w++)
        {
            sum += pSrc[w];
            pDst[w] = sum+pDst_pre[w];
        }
        
        pSrc += step;
        pDst += II_step;
        pDst_pre += II_step;
    }
    
    for (int h=1; h<height; h++)
    {
        sum = 0;
        pDst[-1] = 0;
        int w=0;
        for (; w+4<width; w+=4)
        {
            sum += pSrc[w];
            pDst[w] = sum;
            
            sum += pSrc[w+1];
            pDst[w+1] = sum;
            
            sum += pSrc[w+2];
            pDst[w+2] = sum;
            
            sum += pSrc[w+3];
            pDst[w+3] = sum;
            
            uint32x4_t pre_row = vld1q_u32(pDst_pre+w);
            uint32x4_t cur_row = vld1q_u32(pDst+w);
            cur_row = vaddq_u32(cur_row, pre_row);
            vst1q_u32(pDst+w, cur_row);
        }
        
        for (; w<width; w++)
        {
            sum += pSrc[w];
            pDst[w] = sum+pDst_pre[w];
        }
        
        pSrc += step;
        pDst += II_step;
        pDst_pre += II_step;
    }
    return 0;
}
  #endif
#endif
