SAK
====
This is a base library for other project .     

#### __Add for Ohter Project__    
1. Add `sak` as a submodule in the project    
2. Add following into the project CMakeLists.txt 
```
add_subdirectory(sak)
include_directories(${SAK_INCLUDE_DIR})
```