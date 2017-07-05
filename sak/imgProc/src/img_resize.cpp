#include <math.h>
#include "img_resize.h"


int img_resize(Image_T &src, Image_T &dst)
{
    if (src.format != dst.format)
        return -1;

    switch(src.format)
    {
    case IMG_FMT_GRAY:
        resize_NN_C1(dst.data[0], 
                     dst.pitch[0], dst.width, dst.height,
                     src.data[0], 
                     src.pitch[0], src.width, src.height);
        break;
    case IMG_FMT_BGRBGR:
        resize_NN_C3(dst.data[0], 
                     dst.pitch[0], dst.width, dst.height,
                     src.data[0], 
                     src.pitch[0], src.width, src.height);
        break;
    case IMG_FMT_RGBARGBA:
        resize_NN_C4(dst.data[0], 
                     dst.pitch[0], dst.width, dst.height,
                     src.data[0], 
                     src.pitch[0], src.width, src.height);
        break;
    case IMG_FMT_YYUUVV:
        resize_NN_C1(dst.data[0], 
                     dst.pitch[0], dst.width, dst.height,
                     src.data[0], 
                     src.pitch[0], src.width, src.height);
        resize_NN_C1(dst.data[1], 
                     dst.pitch[1], dst.width, dst.height,
                     src.data[1], 
                     src.pitch[1], src.width, src.height);
        resize_NN_C1(dst.data[2], 
                     dst.pitch[2], dst.width, dst.height,
                     src.data[2], 
                     src.pitch[2], src.width, src.height);
        break;
    default:
        return -1;
    }
    return 0;
}

int resize_NN_C1(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori)
{
    if ((0==img_dst) || (0==img_ori))
        return -1;

    float scale_w = width_ori*1.0f/width_dst;
    float scale_h = height_ori*1.0f/height_dst;

    float src_h=0; 
    for (int h=0; h<height_dst; h++)
    {
        unsigned char *ptr=img_ori+int(src_h+0.5f)*step_ori;
        src_h += scale_h;

        float src_w=0;
        for (int w=0; w<width_dst; w++)
        {
            *(img_dst++) = ptr[int(src_w+0.5f)];
            src_w += scale_w;
        }
    }
    return 0;
}


int resize_NN_C3(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori)
{
    if ((0==img_dst) || (0==img_ori))
        return -1;

    float scale_w = width_ori*1.0f/width_dst;
    float scale_h = height_ori*1.0f/height_dst;

    float src_h=0;
    unsigned char *pDst = img_dst;
    for (int h=0; h<height_dst; h++)
    {
        unsigned char *ptr = img_ori+int(src_h+0.5f)*step_ori;
        src_h += scale_h;
        float src_w = 0;
        for (int w=0; w<width_dst; w++)
        {
            int idx_s = int(src_w+0.5f)*3;
            int idx_d = w*3;
            pDst[idx_d] = ptr[idx_s];
            pDst[idx_d+1] = ptr[idx_s+1];
            pDst[idx_d+2] = ptr[idx_s+2];
            src_w += scale_w;
        }
        pDst += step_dst;
    }
    return 0;
}

int resize_NN_C4(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori)
{
    if ((0==img_dst) || (0==img_ori))
        return -1;

    float scale_w = width_ori*1.0f/width_dst;
    float scale_h = height_ori*1.0f/height_dst;

    float src_h=0;
    unsigned char *pDst = img_dst;
    for (int h=0; h<height_dst; h++)
    {
        unsigned char *ptr = img_ori+int(src_h+0.5f)*step_ori;
        src_h += scale_h;
        float src_w = 0;
        for (int w=0; w<width_dst; w++)
        {
            int idx_s = int(src_w+0.5f)*4;
            int idx_d = w*4;
            pDst[idx_d] = ptr[idx_s];
            pDst[idx_d+1] = ptr[idx_s+1];
            pDst[idx_d+2] = ptr[idx_s+2];
            pDst[idx_d+3] = ptr[idx_s+3];
            src_w += scale_w;
        }
        pDst += step_dst;
    }
    return 0;
}
