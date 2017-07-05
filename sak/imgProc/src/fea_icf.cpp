#include "fea_icf.hpp"

#ifdef ENABLE_NEON
#include "arm_neon.h"

#define SIMD_V4(a) {a, a, a, a}
static uint32x4_t K32_00000001 = SIMD_V4(0x00000001);
static uint32x4_t K32_00000002 = SIMD_V4(0x00000002);
static uint32x4_t K32_00000004 = SIMD_V4(0x00000004);
static uint32x4_t K32_00000008 = SIMD_V4(0x00000008);
static uint32x4_t K32_00000010 = SIMD_V4(0x00000010);
static uint32x4_t K32_00000020 = SIMD_V4(0x00000020);
static uint32x4_t K32_00000040 = SIMD_V4(0x00000040);
static uint32x4_t K32_00000080 = SIMD_V4(0x00000080);

int Fea_MBLBP_4P_Neon(const FeaType &F,
                           const unsigned int *II,
                           const int W, const int H, unsigned int *respond)
{
    const unsigned int *pLine;
    int pitch = F.eh*IntegralImage_Pitch(W);
    
    pLine = IntegralImage_Line(II, W, F.y) + F.x;
    int w1=F.ew;
    int w2=w1<<1;
    int w3=w2+w1;
    uint32x4_t p1[4], p2[4], p3[4], p4[4];
    p1[0] = vld1q_u32(pLine);
    p1[1] = vld1q_u32(pLine+w1);
    p1[2] = vld1q_u32(pLine+w2);
    p1[3] = vld1q_u32(pLine+w3);
    pLine +=  pitch;
    p2[0] = vld1q_u32(pLine);
    p2[1] = vld1q_u32(pLine+w1);
    p2[2] = vld1q_u32(pLine+w2);
    p2[3] = vld1q_u32(pLine+w3);
    pLine +=  pitch;
    p3[0] = vld1q_u32(pLine);
    p3[1] = vld1q_u32(pLine+w1);
    p3[2] = vld1q_u32(pLine+w2);
    p3[3] = vld1q_u32(pLine+w3);
    pLine +=  pitch;
    p4[0] = vld1q_u32(pLine);
    p4[1] = vld1q_u32(pLine+w1);
    p4[2] = vld1q_u32(pLine+w2);
    p4[3] = vld1q_u32(pLine+w3);
    
    //calculate the centre line first
    uint32x4_t val1, val2, val_s, cent;
    uint32x4_t res=vdupq_n_u32(0);
    val1 = vsubq_u32 (p2[0], p3[0]);
    val2 = vsubq_u32 (p2[1], p3[1]);
    val_s = vsubq_u32(val1, val2);
    
    // second line
    // 8
    val1 = vsubq_u32(p2[2], p3[2]);
    cent = vsubq_u32(val2, val1);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000008));
    
    // 16
    val2 = vsubq_u32 (p2[3], p3[3]);
    val_s = vsubq_u32(val1, val2);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000010));
    
    // First line
    //1
    val1 = vsubq_u32(p1[0], p2[0]);
    val2 = vsubq_u32(p1[1], p2[1]);
    val_s = vsubq_u32(val1, val2);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000001));
    //2
    val1 = vsubq_u32(p1[2], p2[2]);
    val_s = vsubq_u32(val2, val1);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000002));
    //4
    val2 = vsubq_u32(p1[3], p2[3]);
    val_s = vsubq_u32(val1, val2);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000004));
    
    // Third line
    //32
    val1 = vsubq_u32(p3[0], p4[0]);
    val2 = vsubq_u32(p3[1], p4[1]);
    val_s = vsubq_u32(val1, val2);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000020));
    //64
    val1 = vsubq_u32(p3[2], p4[2]);
    val_s = vsubq_u32(val2, val1);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000040));
    //128
    val2 = vsubq_u32(p3[3], p4[3]);
    val_s = vsubq_u32(val1, val2);
    res = vorrq_u32(res, vandq_u32(vcgtq_u32(val_s, cent), K32_00000080));
    
    vst1q_u32(respond, res);
    return 0;
}
#endif
