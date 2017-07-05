#ifndef __IMAGE_RESIZE_HPP__
#define __IMAGE_RESIZE_HPP__

#include "image.h"
#ifdef __cplusplus
extern "C" {
#endif

int resize_NN_C1(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori);

int resize_NN_C3(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori);

int resize_NN_C4(unsigned char *img_dst,
                 int step_dst, int width_dst, int height_dst,
                 unsigned char *img_ori,
                 int step_ori, int width_ori, int height_ori);
    
int img_resize(Image_T &src, Image_T &dst);

#ifdef __cplusplus
}
#endif

#endif
