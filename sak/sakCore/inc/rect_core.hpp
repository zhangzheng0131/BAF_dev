#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>
#include "comdef.h"

Rect_T intersect(const Rect_T &r1, const Rect_T &r2);

/*
  [rectList]: 
      Contain the rects to be grouped. And also used to save the result
  [groupThreshold] 
      The min number of each group 
  [eps]
      Use to judge is two rects is similar
  [weights] -- Optional
      Can be an empty vector. If set it will restore the rect's number of each group
  [pConfs] -- Optional
      Contain the conf of each rect. And also restore the conf of grouped rects
 */
void groupRects(std::vector<Rect_T>& rectList, 
                int groupThreshold, double eps, 
                std::vector<int>* weights, 
                std::vector<float>* pConfs);
#endif
