
#最低版本
cmake_minimum_required(VERSION 2.8)
 
#项目信息
project(BSpline)
 
#查找当前目录下所有源文件
#并将文件名保存在DIR_SOURECE
aux_source_directory(. DIR_SOURCE)
 

add_executable(${PROJECT_NAME} ${DIR_SOURCE})

find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
if(GLUT_FOUND)
    include_directories(${GLUT_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${GLUT_LIBRARIES})
else(GLUT_FOUND)
    message(FATAL_ERROR ”GLUT library not found”)
endif(GLUT_FOUND)

if(OPENGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES})
else(OPENGL_FOUND)
    message(FATAL_ERROR ”OpenGL library not found”)
endif(OPENGL_FOUND)

#指定生成目标
