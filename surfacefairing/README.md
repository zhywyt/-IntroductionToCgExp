# 基于余切权值的曲面光顺

## 简介
**项目名: 基于OpenGL 的 余切曲面光顺算法实现。**

**作者: zhywyt**

**日期: 2023-12-17**

**使用了 iGame 的网格数据结构，将 MeshKernel类拿来使用，并自己使用 glut 实现了三角网格的绘制。在`IO.h`中加入了头文件`#include <float.h>`解决宏定义未定义的问题。你可以在`main.cpp` 的`myAlgorithim`函数中实现网格算法，而不需要理会其他的渲染部分。希望能对你带来一些帮助。**

## 编译

### linux

测试环境：
Ubuntu 22.04.3 LTS

- 安装 freeglut
```bash
sudo apt install freeglut3
sudo apt install freeglut3-dev
```
- 安装 Eigen3
```bash
sudo apt install libeigen3-dev
```
#### 2.cmake
```bash
mkdir build
cd build
cmake ..
make
```
### windows
**windows 版本的 Eigen 我暂时不会用，使用 windows 的可以自行配置 Eigen。**


**注：出现任何问题均可联系我，我尽我所能帮你解决编译上的问题。**

**Mainmail: zhywyt@yeah.net**

**HDUmail: zhywyt@hdu.edu.cn**

## 运行
`./demo`
## 用法
- 1 ：切换为`FLAT`光照
- 2 ：切换为`SMOOTH`光照
- a ：相机前进
- A ：相机后退
- GLUT_KEY_LEFT ：物体旋转
- GLUT_KEY_RIGHT ：物体旋转
- GLUT_KEY_UP ：物体旋转
- GLUT_KEY_DOWN ：物体旋转
- j ：执行算法（光顺）