include_directories(./inc) 
aux_source_directory(./src src)

set(target_name "sakImgProc") 
add_library(${target_name} STATIC ${src})
set_property(TARGET ${target_name} PROPERTY POSITION_INDEPENDENT_CODE 1)

