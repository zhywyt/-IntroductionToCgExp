/*
你可以使用左键绘制顶点
你可以使用右键删除顶点
直线的绘制沿着顶点顺序
by zhywyt
最后一次修改
2023.11.26 05:23
*/
#include <iostream>
#include <vector>
#include "DrawLine.hpp"
#include "DrawRound.hpp"
#include "Bspline.hpp"
#include <string>
#include <sstream>
//g++ -o demo main.cpp DrawLine.hpp DrawRound.hpp Point.hpp Bspline.hpp -lGL -lGLU -lglut
using namespace std;
#define m_POINT_SIZE 10
#define m_LINE_SIZE 2
#define WINDOW_HEIGHT 2000
#define WINDOW_WIDTH 1600
typedef enum {
	MOVE,
	CHOOSE,
	NONE
}MOUSEMODE;

vector<myPoint>Vertex;
Circle myCircle;
Ellipt myEllipt;
MOUSEMODE MouseMode = NONE;
GLint ctrlPoint = 0;
BSpline myBSpline;
GLint myBSplineControlIndex = -1,myBSplineControlTemp=-1;

myPoint firstMovePoint{-1,-1},secondMovePoint{-1,-1},temp{-1,-1};
void onDisplay();
void onReshape(GLint w, GLint h);
void onMouse(GLint button, GLint state, GLint x, GLint y);
void onMouseMove(GLint xMouse, GLint yMouse);

void onReshape(GLint w, GLint h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, w, h, 0);
	glMatrixMode(GL_MODELVIEW);
	//发送重绘
	glutPostRedisplay();
}
void onMouse(GLint button, GLint state, GLint x, GLint y){
	printf("Mouse:(%d,%d)\n",x,y);
	if (state == GLUT_DOWN) {
		if (button == GLUT_LEFT_BUTTON) {
			if(MouseMode == CHOOSE){
				myBSplineControlIndex = myBSpline.getPoint(x, y);
				if(myBSplineControlIndex != -1){
					myBSplineControlTemp = myBSplineControlIndex;
					temp = myPoint(myBSpline[myBSplineControlIndex]);
					printf("Control Point: (%d,%d)\n", myBSpline[myBSplineControlIndex].x, myBSpline[myBSplineControlIndex].y);
					}
				else {
					myBSplineControlIndex = myBSplineControlTemp;
				}
				firstMovePoint.x = x;
				firstMovePoint.y = y;
				cout<<"first:"<<firstMovePoint.x<<" "<<firstMovePoint.y<<endl;
			}
			else{
				if(!Vertex.size())
				Vertex.push_back(myPoint(x, y));
				myBSpline.pushControlP(myPoint(x,y));
				ctrlPoint = Vertex.size();
				MouseMode = MOVE;
			}
		}
		else if (button == GLUT_RIGHT_BUTTON) {
			if (MouseMode == MOVE) {
				MouseMode = NONE;
				return;
			}
			auto ibeg = Vertex.begin();
			while (ibeg != Vertex.end()) {
				if (((x - ibeg->x) * (x - ibeg->x)) + ((y - ibeg->y) * (y - ibeg->y)) < 400) {
					Vertex.erase(ibeg);
					break;
				}
				ibeg++;
			}
		}
	}
	if(state == GLUT_UP){
		if(MouseMode == CHOOSE){
			temp = myPoint(myBSpline[myBSplineControlIndex]);
			firstMovePoint .x = -1, firstMovePoint.y = -1;
			secondMovePoint.x = -1, secondMovePoint.y = -1;
		}
	}
	glutPostRedisplay();
}
void onMouseMove(GLint xMouse, GLint yMouse) {
	if (MouseMode == MOVE) {
		Vertex[ctrlPoint-1].x = xMouse, Vertex[ctrlPoint-1].y = yMouse;
		myBSpline[myBSpline.m_n] = Vertex[ctrlPoint-1];
	}
	if(MouseMode == CHOOSE){
		if(myBSplineControlIndex!=-1&&firstMovePoint.x!=-1){
			secondMovePoint.x = xMouse, secondMovePoint.y = yMouse;
			cout<<"first:"<<firstMovePoint.x<<" "<<firstMovePoint.y<<"\n second:"<<secondMovePoint.x<<" "<<secondMovePoint.y<<endl;
			myBSpline[myBSplineControlIndex] = myPoint(secondMovePoint - firstMovePoint + temp);
		}
	}
	glutPostRedisplay();

}
void onMouseActivateMove(GLint xMouse, GLint yMouse){
	if(MouseMode == CHOOSE){
		if(myBSplineControlIndex!=-1&&firstMovePoint.x!=-1){
			secondMovePoint.x = xMouse, secondMovePoint.y = yMouse;
			cout<<"first:"<<firstMovePoint.x<<" "<<firstMovePoint.y<<"\n second:"<<secondMovePoint.x<<" "<<secondMovePoint.y<<endl;
			myBSpline[myBSplineControlIndex] = myPoint(secondMovePoint - firstMovePoint + temp);
		}
	}
	glutPostRedisplay();
}
/*
q&Q 可以退出或者且切换为加点模式
d&D 可以清除所有点
c&C 可以进入选点模式
s&S 可以保存B-Spline文件
*/
void onKeyboard(unsigned char key, int x, int y){
	switch (key) {
	case 'q':
	case 'Q':
		if(MouseMode == CHOOSE){
			cout<<"Move Mode"<<endl;
			MouseMode = MOVE;
			break;
		}
		exit(0);
		break;
	case 'd':
	case 'D':
		if(MouseMode == MOVE||MouseMode==NONE){
			Vertex.clear();
			myBSpline.clear();
			myBSplineControlIndex = -1,myBSplineControlTemp=-1;

			firstMovePoint =myPoint{-1,-1},secondMovePoint=myPoint{-1,-1},temp=myPoint{-1,-1};
		}
		else if(MouseMode == CHOOSE&&myBSplineControlIndex!=-1){
			myBSpline.erase(myBSplineControlIndex);
			myBSplineControlIndex=-1,myBSplineControlTemp=-1;
		}
		break;
	case 'c':
	case 'C':
		//Choose
		cout<<"Choose Mode"<<endl;
		MouseMode = CHOOSE;
		break;
	case 's':
	case 'S':
		//save
		{
			string filename;
			cout<<"input a file name to save the data."<<endl<<":";
			cin>>filename;
			if(filename.find(".txt",0)==string::npos){
				filename += string(".txt");
			}
			if(myBSpline.save(filename)){
				cout<<"save successfully!"<<endl;
				cout<<"at file :"<<filename<<endl;
			}
			else{
				cout<<"save fault!"<<endl;
			}
		}
		//end
		break;
	case 'l':
	case 'L':
		//load
		{
			string filename;
			cout<<"input a file name to load."<<endl<<":";
			cin>>filename;
			if(filename.find(".txt",0)==string::npos){
				filename += string(".txt");
			}
			if(myBSpline.load(filename)){
				cout<<"Load "<<filename<<" successfully!"<<endl;
				myBSpline.showInfo(cout);
			}
			else{
				cout<<"load "<<filename<<" false."<<endl;
			}
		}
		break;
	case 'o':
	case 'O':
		//option
		{
			cout<<"set successfully."<<"\n";
			myBSpline.showInfo(cout);
			
			cout<<"p: set the order of the BSpline"<<endl;
			cout<<"h: set the number of duplicate nodes of the first and last node"<<endl;
			cout<<"input option to set para:"<<endl;
			setbuf(stdin,NULL);
			char opt;
			GLint data;
			cin>>opt;
			cout<<"Input data:"<<endl;
			cin>>data;
			if(myBSpline.setPara(opt,data)){
				cout<<"set successfully."<<"\n";
				myBSpline.showInfo(cout);	
			}
		}
		break;
	case 'h':
	case 'H':
		//help
		cout<<"q Quit\n"
		<<"d Delete\n"
		<<"c Choose\n"
		<<"s Save\n"
		<<"o option\n"
		<<"h help\n"
		<<endl;
		
		break;
	}
}
void onDisplay() {
	glClearColor(224 / 255.0, 237 / 255.0, 253 / 255.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(1.0f, 0, 0);
	auto ibeg = Vertex.begin();
	GLint VertexNum = Vertex.size();
	glPointSize(m_POINT_SIZE);
	glBegin(GL_POINTS);

	while (ibeg != Vertex.end()) {

		glVertex2i(ibeg->x, ibeg->y);
		ibeg++;
	}
	if(firstMovePoint.x!=-1&&secondMovePoint.x!=-1){
		glVertex2i(firstMovePoint.x, firstMovePoint.y);
		glVertex2i(secondMovePoint.x, secondMovePoint.y);
		glPointSize(m_LINE_SIZE);
		BRESENHAM_Line(firstMovePoint,secondMovePoint);
	}
	glEnd();
	if(myBSplineControlIndex!=-1){
		glPointSize(m_LINE_SIZE);
		myCircle.setData(myBSpline[myBSplineControlIndex].x,myBSpline[myBSplineControlIndex].y,m_POINT_SIZE*3);
		myCircle.DrawWithBresenham();
	}

	myBSpline.draw();



	glutSwapBuffers();
	
	// cout << "Once" << endl;
}

GLint main(GLint argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(WINDOW_HEIGHT,WINDOW_WIDTH);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("B-Spline");
	glutReshapeFunc(onReshape);
	glutDisplayFunc(onDisplay);
	glutMouseFunc(onMouse);
	glutPassiveMotionFunc(onMouseMove);
	glutMotionFunc( onMouseActivateMove);
	glutKeyboardFunc(onKeyboard);
	// 设置菜单
	glutMainLoop();
	return 0;
}

