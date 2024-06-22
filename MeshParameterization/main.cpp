#define GLUT_DISABLE_ATEXIT_HACK
#define GL_SILENCE_DEPRECATION
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/freeglut.h>
#include <vector>
#include "IO.h"                                         //iGame框架接口头文件
#include <string>
#include <map>
#include <Eigen/SparseCore>
#include <omp.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
//https://www.licc.tech/article?id=39
using namespace MeshKernel;                             //iGame框架的名字空间
using namespace std;
using namespace Eigen;
GLenum  G_shadingMode;                                  //控制渲染方式
 
float G_fDistance = 3.6f;
float G_fAngle_horizon = 0.0;
float G_fAngle_vertical = 0.0f;
bool grid = false;
 
float G_vLit0Position[4] = { 5.0f, 0.0f, 5.0f, 1.0f };
float G_vLit0Ambient[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit0Diffuse[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit0Specular[4] = { 0.5f, 0.5f, 0.5f, 1.0f };
float G_vMaterialSpecu[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit1Position[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
 
void myinit(void);
void myReshape(GLsizei w, GLsizei h);
void display(void);
void myAlgorithim(GLdouble rate=1e-3);                  //迭代速度
void processSpecialKeys(int key, int x, int y);
void processNormalKeys(unsigned char key,int x,int y);
void drawMyObj();
void showMeshInfo();
void BoundaryRound();
void BoundaryRect();
void LoadMesh();
void global_minimal_surface();

SurfaceMesh myModel;                                    //控制的网格

int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA |GLUT_DEPTH);
    glutInitWindowSize (2000, 2000);
    glutInitWindowPosition (100, 100);
    glutCreateWindow ("OpenGL");

    myinit();

    glutReshapeFunc(myReshape);
    glutSpecialFunc(processSpecialKeys);
    glutKeyboardFunc(processNormalKeys);
    glutDisplayFunc(display);
    glutMainLoop();
    
    return 0;
}
void myinit(void)
{

    glClearColor(0.0, 0.0, 0.0, 0.0);  //背景色
    glShadeModel(GL_FLAT);           
    GLfloat Shininess = 50.0; 
    //材质属性
    glMaterialfv(GL_FRONT, GL_SPECULAR, G_vMaterialSpecu);
    glMaterialfv(GL_FRONT, GL_SHININESS, &Shininess);
    //灯光设置
    glLightfv(GL_LIGHT0, GL_POSITION, G_vLit0Position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, G_vLit0Diffuse);   //散射光属性
    glLightfv(GL_LIGHT0, GL_SPECULAR, G_vLit0Specular);  //镜面反射光
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, G_vLit0Specular);  //环境光参数

    glEnable(GL_LIGHTING);   //开关:使用光
    glEnable(GL_LIGHT0);     //打开0#灯
    glEnable(GL_DEPTH_TEST);

    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);   
    
    //读入网格，这里可以自己 DIY
    IO loader;
    myModel = loader.ReadObjFile("../Model/Nefertiti_face.obj");
    showMeshInfo();
}
void LoadMesh(){
    cout<<"Input the model name:";
    std::string filename;
    cin>>filename;
    filename = "../Model/"+filename+".obj";
    IO loader;
    SurfaceMesh temp = loader.ReadObjFile(filename);
    if(temp.VertexSize()!=0){
        swap(myModel,temp);
    }
}
void showMeshInfo(){
    printf("Succes to load model!\n");
    printf("Face :%ld\n",myModel.FaceSize());
    printf("Vertex :%ld\n",myModel.VertexSize());
    printf("Edge :%ld\n",myModel.EdgeSize());
    printf(myModel.isTriangleMesh()?"Triangle mesh":"Non triangle mesh\n");
}
void myReshape(GLsizei w, GLsizei h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0*(GLfloat)w/(GLfloat)h, 1.0, 30.0);
}
void display(void)
{
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    glClearDepth(0.8f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glShadeModel(G_shadingMode);
    
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    
    glTranslatef(0.0, 0.0, -G_fDistance);
    glRotatef(G_fAngle_horizon, 0.0f, 1.0f, 0.0f);
    glRotatef(G_fAngle_vertical, 1.0f, 0.0f, 0.0f);
    
    //Draw my Obj
    if(grid)glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_TRIANGLES);
    myModel.genAllVerticesNormal();

    auto faces = myModel.allfaces();
    auto vertices = myModel.allvertices();
    int i=0;
    for (auto face : faces) {
        auto v1 = vertices[face.second.vh(0)];
        auto v2 = vertices[face.second.vh(1)];
        auto v3 = vertices[face.second.vh(2)];
        const auto normal = v1.getNormal();
        glNormal3d(normal.x(), normal.y(), normal.z());
        glVertex3f(v1.x(), v1.y(), v1.z());
        const auto normal2 = v2.getNormal();
        glNormal3d(normal2.x(), normal2.y(), normal2.z());
        glVertex3f(v2.x(), v2.y(), v2.z());
        const auto normal3 = v3.getNormal();
        glNormal3d(normal3.x(), normal3.y(), normal3.z());
        glVertex3f(v3.x(), v3.y(), v3.z());
        i++;
    }

    glEnd();
    
    glutSwapBuffers();
}
void processSpecialKeys(int key, int x, int y)
{
    switch(key) {
        case GLUT_KEY_LEFT:
            G_fAngle_horizon -= 5.0f;
            break;
        case GLUT_KEY_RIGHT:
            G_fAngle_horizon += 5.0f;
            break;
        case GLUT_KEY_UP:
            G_fAngle_vertical -= 5.0f;
            break;
        case GLUT_KEY_DOWN:
            G_fAngle_vertical += 5.0f;
            break;
    }
    glutPostRedisplay();
}
void processNormalKeys(unsigned char key,int x,int y)
{
    switch(key) {
        case 'a':  
            G_fDistance -= 0.3f;
            break;
        case 'A':     
            G_fDistance += 0.3f;
            break;
        case '1':   
            G_shadingMode = GL_FLAT;
            break;
        case '2':    
            G_shadingMode = GL_SMOOTH;
            break;
        case 27:    //"esc"
            exit(0);
            break;
        case 'i':
        case 'I':
            LoadMesh();
            break;
        case 'j':
            myAlgorithim();
            showMeshInfo();
            break;
        case 'g':
            global_minimal_surface();
            showMeshInfo();
            break;
        case 'o':
        case 'O':
            BoundaryRound();
            break;
        case 'R':
        case 'r':
            BoundaryRect();
            break;
        case 'l':
        case 'L':
            grid=!grid;
    }
    glutPostRedisplay();
}
void myAlgorithim(GLdouble rate){
    //获取所有的顶点
    auto vertices = myModel.allvertices();
    // 创建一个新的顶点集合来存储更新后的顶点位置
    std::map<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex> newVertices;
    for(auto &v : vertices){
        //如果不是边界点，则需要更新该点
        if(!myModel.isOnBoundary(v.first)){
            //点的结构
            //      v1       //
		    //    //||\\     //
            //   v2 || v3    //
            //    \\||//     //
            //      v0       //
            //获得 v0 相连的所有边
            auto neb = myModel.NeighborEh(v.first);
            MeshKernel::iGameVertex lamplacian;
            double area = 0;
            //遍历相连的边
            for(auto eh:neb){
                const auto v0h = v.first;
                //通过一个点和一个边，找到另一个点
                const auto v1h = myModel.NeighborVhFromEdge(v0h,eh);
                MeshKernel::iGameFaceHandle f1,f2;
                MeshKernel::iGameVertexHandle v2h,v3h;
                //找到 v0 v1 这条边的相邻面
                const auto nebface = myModel.NeighborFh(eh);
                int i=0;
                //取出两个面（流形只会有两个相邻面）
                for(auto fh:nebface){
                    if(i==0){
                        f1 = fh;
                        i++;
                    }
                    else{
                        f2 = fh;
                    }
                }
                //找到这两个面上另外的两个点 v2 v3
                for(i=0;i<3;i++){
                    if(myModel.faces(f1).vh(i)!=v0h&&myModel.faces(f1).vh(i)!=v1h)
                        v2h = myModel.faces(f1).vh(i); 
                    if(myModel.faces(f2).vh(i)!=v0h&&myModel.faces(f2).vh(i)!=v1h)
                        v3h = myModel.faces(f2).vh(i);
                }
                //复制这四个点的信息
                const auto v0 = myModel.vertices(v0h);
                const auto v1 = myModel.vertices(v1h);
                const auto v2 = myModel.vertices(v2h);
                const auto v3 = myModel.vertices(v3h);
                //计算需要使用的向量
                const auto dp = v0-v1;
                const auto v21 = v2 - v1;
                const auto v31 = v3 - v1;
                const auto v20 = v2 - v0;
                const auto v30 = v3 - v0;
                // % 操作符是叉乘，计算两个面的面积
                area += (v21%v20).norm()/2+(v31%v30).norm()/2;
                // 使用 cos 比上 sin 计算 cot 
                double cot2 = (v21*v20)/(v21%v20).norm();
                double cot1 = (v31*v30)/(v31%v30).norm();
                //更新 lamplacin 
                lamplacian += dp*(-cot1-cot2);
            }
            //本来是除以 2*area，但是因为每个面计算了两遍，所以不需要处理。
            lamplacian/=area;
            lamplacian*=rate;
            // 在新的顶点集合中存储更新后的顶点位置，而不是直接更新原始模型
            newVertices[v.first] = myModel.vertices(v.first) + lamplacian;
        }
        else {
            // 对于边界顶点，我们不进行更新
            newVertices[v.first] = myModel.vertices(v.first);
        }
    }
    // 现在我们更新原始模型的顶点位置
    for(auto &v : newVertices){
        myModel.vertices(v.first) = v.second;
    }
}
//把边界点映射到圆上
void BoundaryRound(){
    iGameVertex O(0,0,0);
    int num=0;
    auto vertices = myModel.allvertices();
    for(auto&v:vertices){
        if(myModel.isOnBoundary(v.first))O+=v.second,num++;
    }
    std::map<iGameVertexHandle,iGameVertex>newBuondary;
    if(num)
        O/=num;
    else {
        cout<<"There is not boundary!"<<endl;
        return;
    }
    for(auto&v:vertices){
        if(myModel.isOnBoundary(v.first)){
            iGameVertex temp=v.second-O;
            temp.setZ(0);
            temp = temp.normalize();
            newBuondary[v.first]=iGameVertex(temp.x(),temp.y(),0);
        }
    }
    for(auto&v:newBuondary){
        myModel.vertices(v.first)=v.second;
    }
}
void BoundaryRect(){
    iGameVertex Rec(0,0,0);
    int num=0;
    auto vertices = myModel.allvertices();
    for(auto&v:vertices){
        if(myModel.isOnBoundary(v.first))Rec+=v.second,num++;
    }
    std::map<iGameVertexHandle,iGameVertex>newBuondary;
    if(num)
        Rec/=num;
    else {
        cout<<"There is not boundary!"<<endl;
        return;
    }
    for(auto&v:vertices){
        if(myModel.isOnBoundary(v.first)){
            iGameVertex temp=v.second - Rec;
            temp.setZ(0);
            temp = temp.normalize();
            if(abs(temp.x())>abs(temp.y())){
                if(temp.x()>0){
                    temp.setY(temp.y()/temp.x());
                    temp.setX(1);
                }
                else{
                    temp.setY(temp.y()/-temp.x());
                    temp.setX(-1);
                }
            }
            else{
                if(temp.y()>0){
                    temp.setX(temp.x()/temp.y());
                    temp.setY(1);
                }
                else{
                    temp.setX(temp.x()/-temp.y());
                    temp.setY(-1);
                }
            }
            newBuondary[v.first]=iGameVertex(temp.x(),temp.y(),0);
        }
    }
    for(auto&v:newBuondary){
        myModel.vertices(v.first)=v.second;
    }
}

void global_minimal_surface(){
    //数据准备
    Eigen::SparseMatrix<float>A(myModel.VertexSize(),myModel.VertexSize());
    Eigen::SparseLU<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> >solver;
    vector<Eigen::Triplet<float> >tv;
    Eigen::VectorXf b_x(myModel.VertexSize());
    Eigen::VectorXf b_y(myModel.VertexSize());
    Eigen::VectorXf b_z(myModel.VertexSize());
    //填充L矩阵
    auto vertices = myModel.allvertices();
    for(auto v:vertices){
        if(myModel.isOnBoundary(v.first)){
            Eigen::Triplet<float> temp(v.first,v.first,1.0f);
            tv.push_back(temp);
            b_x[v.first] = v.second.x();
            b_y[v.first] = v.second.y();
            b_z[v.first] = v.second.z();
        }
        else
        {
            // 点的结构
            //       v1       //
            //     //||\\     //
            //    v2 || v3    //
            //     \\||//     //
            //       v0       //
            // 获得 v0 相连的所有边
            auto neb = myModel.NeighborEh(v.first);
            float cont = 0;
            // 遍历相连的边
            for (auto eh : neb)
            {
                const auto v0h = v.first;
                // 通过一个点和一个边，找到另一个点
                const auto v1h = myModel.NeighborVhFromEdge(v0h, eh);
                MeshKernel::iGameFaceHandle f1, f2;
                MeshKernel::iGameVertexHandle v2h, v3h;
                // 找到 v0 v1 这条边的相邻面
                const auto nebface = myModel.NeighborFh(eh);
                int i = 0;
                // 取出两个面（流形只会有两个相邻面）
                for (auto fh : nebface)
                {
                    if (i == 0)
                    {
                        f1 = fh;
                        i++;
                    }
                    else
                    {
                        f2 = fh;
                    }
                }
                // 找到这两个面上另外的两个点 v2 v3
                for (i = 0; i < 3; i++)
                {
                    if (myModel.faces(f1).vh(i) != v0h && myModel.faces(f1).vh(i) != v1h)
                        v2h = myModel.faces(f1).vh(i);
                    if (myModel.faces(f2).vh(i) != v0h && myModel.faces(f2).vh(i) != v1h)
                        v3h = myModel.faces(f2).vh(i);
                }
                // 复制这四个点的信息
                const auto v0 = myModel.vertices(v0h);
                const auto v1 = myModel.vertices(v1h);
                const auto v2 = myModel.vertices(v2h);
                const auto v3 = myModel.vertices(v3h);
                // 计算需要使用的向量
                const auto v21 = v2 - v1;
                const auto v31 = v3 - v1;
                const auto v20 = v2 - v0;
                const auto v30 = v3 - v0;
                // 使用 cos 比上 sin 计算 cot
                double cot2 = (v21 * v20) / (v21 % v20).norm();
                double cot1 = (v31 * v30) / (v31 % v30).norm();
                float w = cot1+cot2;
                Eigen::Triplet<float> temp(v0h,v1h,w);
                tv.push_back(temp);
                cont+=w;
            }
            Eigen::Triplet<float>temp(v.first,v.first,-cont);
            tv.push_back(temp);
            b_x[v.first]=0.0f;
            b_y[v.first]=0.0f;
            b_z[v.first]=0.0f;
        }
    }
    //Solve
    A.setFromTriplets(tv.begin(),tv.end());
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXf x = solver.solve(b_x);
    Eigen::VectorXf y = solver.solve(b_y);
    Eigen::VectorXf z = solver.solve(b_z);
    for(auto v:vertices){
        if(myModel.isOnBoundary(v.first))continue;
        myModel.vertices(v.first)=iGameVertex(x[v.first],y[v.first],z[v.first]);
    }
}
