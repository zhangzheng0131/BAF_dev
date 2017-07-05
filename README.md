Learning Background-Aware Correlation Filters for Visual Tracking
====

## Build
```
mkdir build; cd build; cmake ..; make
```

## Supported Algorithms

* KCF "High-Speed Tracking with Kernelized Correlation Filters"
* Staple "Staple: Complementary Learners for Real-Time Tracking"

## Todo List

* Apply gaussian kernel into staple
* Apply LUV color feature into staple


## introduction
Bacf implementation of c++ based on staple. Add function ECF according to Major changes happen at function : trainTransCF.
other function in file stapletracker.cpp and demos.cpp (eg: demo_img_list.cpp) are all changed more or less to fit for
implementation of bacf.

for viewing the farmework of bacf. you can see my annotation which mainly correspond to matlab implementation of bacf.you can compare differences between them to understand each other.


## tips for running on OTB100
* point  path to ./BACFcpp/benchmark/ in main_running.m
I implement like this:
                //
                switch t.name
                    case {'VR','TM','RS','PD','MS'}
*                    case {'BACFcpp'}
                        cd(['./trackers/' t.name '/benchmark/']);
                        addpath(genpath('./'))
                    otherwise
                        cd(['./trackers/' t.name]);
                        addpath(genpath('./'))
                end

                eval(funcName);

                switch t.name
                    case {'VR','TM','RS','PD','MS'}
*                    case {'BACFcpp'}
                        rmpath(genpath('./'))
                        cd('../../../');
                    otherwise
                        rmpath(genpath('./'))
                        cd('../../');
                end
                //
