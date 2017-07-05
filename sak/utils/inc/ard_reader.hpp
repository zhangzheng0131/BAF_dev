#ifndef __ARD_READER_HPP__
#define __ARD_READER_HPP__

#include <string>
#include <vector>
#include <stdio.h>
#include "image.h"

class ARDReader
{
public:
    ARDReader(std::string path)
    {
        m_ardPath = path;
        m_ardFile = 0;
        m_format  = "";
        std::string fmt[] = {"YYUUVV",
                             "BGRBGR",
                             "GRAY"}; 
        m_validFMT.assign(fmt, 
                          fmt+sizeof(fmt)/sizeof(fmt[0]));
    }
    ~ARDReader()
    {
        uninit();
    }

public:
    int init();
    int nextFrame(Image_T &img);    
    inline int channel(){return m_ch;}
    inline int width(){return m_width;}
    inline int height(){return m_height;}
    
private:
    void uninit();
    bool isValidFMT();
    
private:
    FILE *m_ardFile;
    std::string m_ardPath;
    int m_width;
    int m_height;
    int m_ch;
    std::string m_format;
    std::vector<std::string> m_validFMT;
};




#endif
