# 曲面参数化

## 简介
**项目名: 基于OpenGL 的 曲面参数化实现。**

**作者: zhywyt**

**日期: 2024-01-07**


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

## 文件说明

### main.cpp
```cpp
void myAlgorithim(GLdouble rate=1e-3);                  //曲面光顺
void processSpecialKeys(int key, int x, int y);         //键盘输入
void processNormalKeys(unsigned char key,int x,int y);  //键盘输入
void BoundaryRound();                                   //圆形边界映射  
void BoundaryRect();                                    //矩形边界映射
void LoadMesh();                                        //导入网格
void global_minimal_surface();                          //全局曲面光顺
```
**在对边界点进行映射之后，再执行全局曲面光顺就可以达到参数化的效果，但其实会存在三角形反向的情况。**

## 运行
`./demo`
## 用法
- 1 : 切换为`FLAT`光照
- 2 : 切换为`SMOOTH`光照
- l : 切换为相框模式
- a : 相机前进
- A : 相机后退
- GLUT_KEY_LEFT : 物体旋转
- GLUT_KEY_RIGHT : 物体旋转
- GLUT_KEY_UP : 物体旋转
- GLUT_KEY_DOWN : 物体旋转
- j : 执行算法（光顺）
- g : 执行算法 (参数化)
- i : 导入obj文件,通过控制台输入模型名字实现
- o : 将边界点映射到圆形区域
- r : 将边界点映射到矩形区域

