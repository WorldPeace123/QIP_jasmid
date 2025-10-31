#
# 文件名: examples/LinAdvSL/source/CMakeLists.txt
# 软件包: JASMIN cmake system
# 版权  : (c) 2004-2017 北京应用物理与计算数学研究所
# 版权  : (c) 2013-2017 中物院高性能数值模拟软件中心
# 版本号: $Revision$
# 修改  : $Date$
# 描述  : 例子程序的CMakeLists.txt
#
# 使用示例:
#
# 方法一: 使用cmake和make
# $ cd build/
# $ cmake -D JASMIN_NDIM=2|3 ..
# $ make main2d|main3d
#
# 方法二: 使用make（和3.0相同）
# $ cd source/
# $ make main-2d
# $ or
# $ make main-3d

cmake_minimum_required(VERSION 3.3)
#find_package(MPI REQUIRED)
project(JASMIN_example CXX C Fortran)

if(NOT JASMIN_ROOT)
  set(JASMIN_ROOT /work1/luteng/soft/jusmid)
endif()

if(NOT JASMIN_NDIM)
  set(JASMIN_NDIM 3)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE )
endif()

list(APPEND CMAKE_MODULE_PATH ${JASMIN_ROOT}/share/cmake/Modules/scns)
find_package(JASMIN REQUIRED)
add_definitions(${JASMIN_DEFINITIONS})
include_directories(${JASMIN_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR})

# 源文件
set(SOURCES
    main_QIP.C
    QIP.C
    QIPLevelIntegrator.C
)

# 头文件
set(HEADERS
    QIP.h
    QIPLevelIntegrator.h
)
set(CMAKE_BUILD_TYPE Debug)
# 创建可执行文件
add_executable(qip_jasmin ${SOURCES} ${HEADERS})

# 链接库
target_link_libraries(qip_jasmin ${JASMIN_LIBRARIES})

# 编译选项
target_compile_options(qip_jasmin PRIVATE ${JASMIN_COMPILE_FLAGS})

# 安装
install(TARGETS qip_jasmin DESTINATION bin)

