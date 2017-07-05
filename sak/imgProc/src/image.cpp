#include "image.h"
#include <string.h>
int allocImage(Image_T &img, int width, int height, 
               IMG_FMT format)
{
    int res = 0;
    memset(&img, 0, sizeof(Image_T));
    img.width = width;
    img.height = height;
    img.format = format;
    switch(format)
    {
    case IMG_FMT_GRAY:
        img.nPlane=1;
        img.data[0] = new unsigned char[width*height];
        img.pitch[0] = width;
        break;
    case IMG_FMT_BGRBGR:
        img.nPlane=1;
        img.data[0] = new unsigned char[width*height*3];
        img.pitch[0] = width*3;
        break;
    case IMG_FMT_YYUUVV:
        img.nPlane=3;
        img.data[0] = new unsigned char[width*height*3];
        img.pitch[0] = img.pitch[1] = img.pitch[2] = width;
        img.data[1] = img.data[0] + width*height;
        img.data[2] = img.data[1] + width*height;
        break;
    case IMG_FMT_UNKNOWEN:
        img.data[0] = new unsigned char[width*height*4];
        break;
    default:
        res = -1;       
    }
    return res;
}

void deallocImage(Image_T &img)
{
    if (0 != img.data[0])
        delete []img.data[0];
    memset(&img, 0, sizeof(Image_T));
    return ;
}

