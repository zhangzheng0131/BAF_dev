
### Set the SAK_INCLUDE_DIR for parent scope
set(SAK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/common
                    ${CMAKE_CURRENT_SOURCE_DIR}/sakCore/inc
                    ${CMAKE_CURRENT_SOURCE_DIR}/utils/inc
                    ${CMAKE_CURRENT_SOURCE_DIR}/imgProc/inc
                    ${CMAKE_CURRENT_SOURCE_DIR}/jsmn/inc
                    PARENT_SCOPE)

### Add include directories for sak
include_directories(./common)
include_directories(./sakCore/inc)
include_directories(./imgProc/inc)
include_directories(./utils/inc)
include_directories(./jsmn/inc)

### Build the sub libraries
add_subdirectory(./utils)
add_subdirectory(./imgProc)
add_subdirectory(./sakCore)
add_subdirectory(./jsmn)


install(FILES ./common/comdef.h 
        DESTINATION inc)
