#include <stdio.h>
#include <string.h>
#include "log.h"

static int optind = 1;
static char *optarg;

int getpos_t()
{
    return optind;
}
char* getarg_t()
{
    return optarg;
}

int getopt_t(int argc, char *argv[], char *opts)
{
    static int sp = 1;
    int c;
    char *cp;

    if (sp == 1) {
        if (optind >= argc ||
            argv[optind][0] != '-' || argv[optind][1] == '\0')
            return EOF;
        else if (!strcmp(argv[optind], "--")) {
            optind++;
            return EOF;
        }
    }
    c = argv[optind][sp];
    if (c == ':' || !(cp = strchr(opts, c))) {
        fprintf(stderr, ": illegal option -- %c\n", c);
        if (argv[optind][++sp] == '\0') {
            optind++;
            sp = 1;
        }
        return '?';
    }
    if (*++cp == ':') {
        if (argv[optind][sp+1] != '\0')
            optarg = &argv[optind++][sp+1];
        else if(++optind >= argc) {
            fprintf(stderr, ": option requires an argument -- %c\n", c);
            sp = 1;
            return '?';
        } else
            optarg = argv[optind++];
        sp = 1;
    } else {
        if (argv[optind][++sp] == '\0') {
            sp = 1;
            optind++;
        }
        optarg = NULL;
    }

    return c;
}

#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <sys/time.h>
#include <time.h>
double timeStamp()
{
    double usec = 0;
    struct timeval time;
    gettimeofday(&time, 0);

    usec = time.tv_sec*1000000 + time.tv_usec;
    return  usec;
}

int isExpired(int EXPIRATION_YEAR, int EXPIRATION_MONTH)
{
    struct timeval tv;  
    struct tm *tm;    
    gettimeofday(&tv, NULL);  
    tm = localtime(&tv.tv_sec); 

    int year  = tm->tm_year + 1900;
    int month = tm->tm_mon + 1;
    //int day   = tm->tm_mday;
    
    if (year > EXPIRATION_YEAR)
        return 1;
    else if (year == EXPIRATION_YEAR && 
             month > EXPIRATION_MONTH)
        return 1;
    return 0;
}
#elif defined(WIN32) || defined(WIN64)
#include <Windows.h>
double timeStamp()
{
    double usec=0;
    LARGE_INTEGER freq, t;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t);
    usec = t.QuadPart*1000000.0 / freq.QuadPart;
    return  usec;
}

int isExpired(int EXPIRATION_YEAR, int EXPIRATION_MONTH)
{
    SYSTEMTIME st;
    GetLocalTime(&st);

    if (st.wYear > EXPIRATION_YEAR)
        return 1;
    else if (st.wYear == EXPIRATION_YEAR && 
             st.wMonth > EXPIRATION_MONTH)
        return 1;
    return 0;   
}
#else
//#error ("Unsupported platform")
#endif

char* strstrip(char *s)
{
    char *p = s;
    int l = strlen(p);
    if (0 == l)
        return s;
    
    while(' '==p[l-1] || '\n'==p[l-1]) 
        p[--l] = '\0';
    while(*p && (' '==*p || '\n'==*p))
    {
        ++p;
        --l;
    }
    memmove(s, p, l+1);
    return s;
}
