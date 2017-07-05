#include <string.h>
#include <math.h>
#include <stdio.h>
#include "comdef.h"
#include "img_resize.h"
#include "img_convert.h"
#include "cimage.hpp"

CImage_T::CImage_T(int maxSide, IMG_FMT format)
{
    m_maxSide = maxSide;
    m_format  = format;
    m_scale = 0.f;
    memset(&m_img, 0, sizeof(Image_T));
    memset(&m_tmp, 0, sizeof(Image_T));
}

CImage_T::~CImage_T()
{
    uninit();
}

int CImage_T::setImage(Image_T &img, float &rate)
{
    m_scale = MAX_T(1.f,
                   MAX_T(img.width*1.f/m_maxSide,
                          img.height*1.f/m_maxSide));
    
    int width  = round(img.width/m_scale);
    int height = round(img.height/m_scale);
    rate=m_scale;
    //int width=img.width;
   //int height=img.height;

    if (0 != resetImage(m_img, width, height, m_format))
        return -1;

    if (0 != resetImage(m_tmp, width, height, img.format))
       return -1;

    if (img.format==m_img.format &&
        width==img.width &&
        height==img.height)
   {
        if (0 != copyImage(img, m_img, m_format)) {
            printf("img is the same as m_img\n");
            return -1;
        }
    }
   else if (img.format == m_img.format)
    {
        if (0 != img_resize(img, m_img)) {
            printf("It seems img should be resize if too big \n");
            return -1;
        }
    }
    else if (width==img.width && height==img.height)
    {
        if (0 != img_convert(img, m_img))
            return -1;
    }
    else
    {
        if (0 != img_resize(img, m_tmp))
            return -1;
        if (0 != img_convert(m_tmp, m_img))
            return -1;
    }

  //  printf("where is the point of question ??\n");
//    copyImage(img, m_img, m_format);
    return 0;
}

int CImage_T::init()
{
    if (0 != allocImage(m_img, m_maxSide, 
                        m_maxSide, m_format))
        return -1;
    if (0 != allocImage(m_tmp, m_maxSide, 
                        m_maxSide, IMG_FMT_UNKNOWEN))
        return -1;
    return 0;
}

void CImage_T::uninit()
{
    deallocImage(m_img);
    deallocImage(m_tmp);
}


int CImage_T::copyImage(Image_T &src, Image_T &dst, 
                       IMG_FMT format)
{
    if (src.width!=dst.width || src.height!=dst.height)
        return -1;
    if (src.format!=format || dst.format!=format)
        return -1;

    switch(format)
    {
    case IMG_FMT_GRAY:
        memcpy(dst.data[0], src.data[0],
               src.width*src.height);
        break;
    case IMG_FMT_BGRBGR:
        memcpy(dst.data[0], src.data[0], 
               src.width*src.height*3);
        break;
    case IMG_FMT_YYUUVV:
        memcpy(dst.data[0], src.data[0], 
               src.width*src.height*3);
        break;
    case IMG_FMT_RGBARGBA:
        memcpy(dst.data[0], src.data[0], 
               src.width*src.height*4);
        break;
    default:
        return -1;       
    }

    return 0;
}

int CImage_T::resetImage(Image_T &img, 
                         int width, int height, 
                         IMG_FMT format)
{
    int res = 0;

    img.width = width;
    img.height = height;
    switch(format)
    {
    case IMG_FMT_GRAY:
        img.nPlane=1;
        img.pitch[0] = width;
        img.format = IMG_FMT_GRAY;
        break;
    case IMG_FMT_BGRBGR:
        img.nPlane=1;
        img.pitch[0] = width*3;
        img.format = IMG_FMT_BGRBGR;
        break;
    case IMG_FMT_RGBARGBA:
        img.nPlane=1;
        img.pitch[0] = width*4;
        img.format = IMG_FMT_RGBARGBA;
        break;
    case IMG_FMT_YYUUVV:
        img.nPlane=3;
        img.pitch[0] = img.pitch[1] = img.pitch[2] = width;
        img.data[1] = img.data[0] + width*height;
        img.data[2] = img.data[1] + width*height;
        img.format = IMG_FMT_YYUUVV;
        break;
    default:
        res = -1;       
    }
    return res;
}
