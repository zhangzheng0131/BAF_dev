#ifndef __CIMAGE_HPP__
#define __CIMAGE_HPP__

#include "image.h"
class CImage_T
{
public:
    CImage_T(int maxSide, IMG_FMT format);
    virtual ~CImage_T();
    
public:
    int setImage(Image_T &img ,float &rate);

protected:
    int init();

private:
    int resetImage(Image_T &img, 
                   int width, int height, IMG_FMT format);
    int copyImage(Image_T &src, Image_T &dst, IMG_FMT format);
    void uninit();

private:
    Image_T  m_tmp;

protected:
    float    m_scale;
    Image_T  m_img;    
    int      m_maxSide;
    IMG_FMT  m_format;
};



#endif
