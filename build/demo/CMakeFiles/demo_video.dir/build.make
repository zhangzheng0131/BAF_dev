# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build

# Include any dependencies generated for this target.
include demo/CMakeFiles/demo_video.dir/depend.make

# Include the progress variables for this target.
include demo/CMakeFiles/demo_video.dir/progress.make

# Include the compile flags for this target's objects.
include demo/CMakeFiles/demo_video.dir/flags.make

demo/CMakeFiles/demo_video.dir/demo_video.cpp.o: demo/CMakeFiles/demo_video.dir/flags.make
demo/CMakeFiles/demo_video.dir/demo_video.cpp.o: ../demo/demo_video.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object demo/CMakeFiles/demo_video.dir/demo_video.cpp.o"
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo && /bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/demo_video.dir/demo_video.cpp.o -c /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/demo/demo_video.cpp

demo/CMakeFiles/demo_video.dir/demo_video.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/demo_video.dir/demo_video.cpp.i"
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo && /bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/demo/demo_video.cpp > CMakeFiles/demo_video.dir/demo_video.cpp.i

demo/CMakeFiles/demo_video.dir/demo_video.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/demo_video.dir/demo_video.cpp.s"
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo && /bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/demo/demo_video.cpp -o CMakeFiles/demo_video.dir/demo_video.cpp.s

demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.requires:
.PHONY : demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.requires

demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.provides: demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.requires
	$(MAKE) -f demo/CMakeFiles/demo_video.dir/build.make demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.provides.build
.PHONY : demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.provides

demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.provides.build: demo/CMakeFiles/demo_video.dir/demo_video.cpp.o

# Object files for target demo_video
demo_video_OBJECTS = \
"CMakeFiles/demo_video.dir/demo_video.cpp.o"

# External object files for target demo_video
demo_video_EXTERNAL_OBJECTS =

bin/demo_video: demo/CMakeFiles/demo_video.dir/demo_video.cpp.o
bin/demo_video: demo/CMakeFiles/demo_video.dir/build.make
bin/demo_video: bin/libsakUtils.a
bin/demo_video: bin/libobject_tracker.a
bin/demo_video: /home/zhangzheng/software/lib/libopencv_videostab.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_video.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_ts.a
bin/demo_video: /home/zhangzheng/software/lib/libopencv_superres.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_stitching.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_photo.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_objdetect.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_nonfree.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_ml.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_legacy.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_imgproc.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_highgui.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_gpu.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_flann.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_features2d.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_core.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_contrib.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_calib3d.so.2.4.11
bin/demo_video: bin/libsakUtils.a
bin/demo_video: /home/zhangzheng/software/lib/libopencv_nonfree.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_gpu.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_photo.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_objdetect.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_legacy.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_video.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_ml.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_calib3d.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_features2d.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_highgui.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_imgproc.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_flann.so.2.4.11
bin/demo_video: /home/zhangzheng/software/lib/libopencv_core.so.2.4.11
bin/demo_video: bin/libsakImgProc.a
bin/demo_video: demo/CMakeFiles/demo_video.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/demo_video"
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo_video.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demo/CMakeFiles/demo_video.dir/build: bin/demo_video
.PHONY : demo/CMakeFiles/demo_video.dir/build

demo/CMakeFiles/demo_video.dir/requires: demo/CMakeFiles/demo_video.dir/demo_video.cpp.o.requires
.PHONY : demo/CMakeFiles/demo_video.dir/requires

demo/CMakeFiles/demo_video.dir/clean:
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo && $(CMAKE_COMMAND) -P CMakeFiles/demo_video.dir/cmake_clean.cmake
.PHONY : demo/CMakeFiles/demo_video.dir/clean

demo/CMakeFiles/demo_video.dir/depend:
	cd /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/demo /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo /home/zhangzheng/workplace/master/nntrackerpapers/matlab/tracker_benchmark_v1.1/tracker_benchmark_v1.1/trackers/BACFcpp/build/demo/CMakeFiles/demo_video.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demo/CMakeFiles/demo_video.dir/depend

