#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ard_reader.hpp"
#include "comdef.h"
#include "log.h"

bool ARDReader::isValidFMT()
{
    std::vector<std::string>::iterator it;
    it = m_validFMT.begin();
    for(; it<m_validFMT.end(); it++)
    {
        if (0 == m_format.compare(*it))
            return true;
    }
    return false;
}

int ARDReader::init()
{
    uninit();
    m_ardFile = fopen(m_ardPath.c_str(), "rb");
    if (0==m_ardFile)
        return -1;
    
    bool isValid = false;
    char flag[1024] = {0};
    fgets(flag, 1024, m_ardFile);
    flag[strlen(flag)-1] = 0;    
    if (0 != strcmp(flag, "<HEADER>"))
        return -1;
    
    while (!feof(m_ardFile))
    {
        memset(flag, 0, 1024);        
        fgets(flag, 1024, m_ardFile);
        flag[strlen(flag)-1] = 0;   
        
        if (0 == strcmp(flag, "<!HEADER>"))
        {
            isValid = true;
            break;
        }
        
        char *val = strstrip(strtok(flag, ":"));
        if (0==strcmp(val, "WIDTH"))
            m_width = atoi(strstrip(strtok(0, ":")));
        else if (0==strcmp(val, "HEIGHT"))
            m_height = atoi(strstrip(strtok(0, ":")));
        else if (0==strcmp(val, "FORMAT"))
            m_format = strstrip(strtok(0, ":"));
    }
    
    if (false == isValid)
        return -1;

    if (false == isValidFMT())
        return -1;

    // Set the channel value
    if (0 == m_format.compare("GRAY"))
        m_ch = 1;
    else if (0 == m_format.compare("YYUUVV"))
        m_ch = 3;
    else if (0 == m_format.compare("BGRBGR"))
        m_ch = 3;
    else
        m_ch = 0;
    return 0;
}

void ARDReader::uninit()
{
    if (0 != m_ardFile)
    {
        fclose(m_ardFile);
        m_ardFile = 0;
    }
}

int ARDReader::nextFrame(Image_T &img)
{
    if (feof(m_ardFile))
        return -1;
    fread(img.data[0], m_width*m_height*m_ch, 1, m_ardFile);
    img.width = m_width;
    img.height = m_height;    
    
    if (0 == m_format.compare("GRAY"))
    {
        img.pitch[0] = m_width;
        img.format = IMG_FMT_GRAY;     
    }
    else if (0 == m_format.compare("YYUUVV"))
    {
        img.pitch[0] = m_width;
        img.pitch[1] = m_width;
        img.pitch[2] = m_width;
        img.data[1] = img.data[0]+m_width*m_height;
        img.data[2] = img.data[1]+m_width*m_height;
        img.format = IMG_FMT_YYUUVV;     
    }
    else if (0 == m_format.compare("BGRBGR"))
    {
        img.pitch[0] = m_width*m_ch;
        img.format = IMG_FMT_BGRBGR;     
    }
    return 0;
}
