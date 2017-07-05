#ifndef __COMDEF_H__
#define __COMDEF_H__

#include <stdint.h>

typedef struct _tagRect
{
    int x;
    int y;
    int w;
    int h;
}Rect_T;

typedef enum 
{
    IMG_FMT_ERROR    = -1,
    IMG_FMT_UNKNOWEN = 0,
    IMG_FMT_GRAY     = 1,
    IMG_FMT_BGRBGR   = 2,
    IMG_FMT_RRGGBB   = 3,
    IMG_FMT_YYUUVV   = 4,
    IMG_FMT_RGBARGBA   = 5
}IMG_FMT;

typedef struct __tagImage
{
    unsigned char *data[4];
    unsigned int  pitch[4];
    int width;
    int height;
    int nPlane;
    IMG_FMT format;
}Image_T;

#define PI_T 3.14159265f

#define MAX_T(a,b)  ((a)>(b)?(a):(b))
#define MIN_T(a,b)  ((a)<(b)?(a):(b))

#endif
