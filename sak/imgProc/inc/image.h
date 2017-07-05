#ifndef __IMAGE_H__
#define __IMAGE_H__

#include "comdef.h"
#ifdef __cplusplus
extern "C" {
#endif

int allocImage(Image_T &img,
               int width, int height, IMG_FMT format);

void deallocImage(Image_T &img);


#ifdef __cplusplus
}
#endif

#endif
