/*
这是绘制圆形的算法
*/


#ifndef DRAWROUND_HPP
#define DRAWROUND_HPP
#include <cmath>
#include <GL/freeglut.h>

#include "Point.hpp"

//圆
class Circle {
public:
	Circle(GLint x=0, GLint y=0, GLint r=0):Rx(x),Ry(y),R(r) {}
	Circle(GLint rx, GLint ry, GLint endx, GLint endy) :Rx(rx),Ry(ry){
		R = static_cast<GLint>(sqrt((endx - rx) * (endx - rx) + (endy - ry) * (endy - ry)));
	}
	void setData(GLint rx, GLint ry, GLint endx, GLint endy) {
		Rx = rx, Ry = ry;
		R = static_cast<GLint>(sqrt((endx - rx) * (endx - rx) + (endy - ry) * (endy - ry)));
	}
	void setData(GLint x , GLint y , GLint r ) {
		Rx = x, Ry = y, R = r;
	}
	void DrawWithMid() {
		Mid_Circle(Rx, Ry, R);
	}
	void DrawWithBresenham() {
		BRESENHAM_Circle(Rx,Ry,R);
	}
private:
	GLint Rx, Ry, R;
	
	/*
	中点画圆法
	*/
	/// <summary>
	/// 中点画圆算法
	/// </summary>
	/// <param name="Rx">圆心</param>
	/// <param name="Ry"></param>
	/// <param name="R">半径</param>
	void Mid_Circle(GLint Rx, GLint Ry, GLint R) {
		glBegin(GL_POINTS);
		GLint x, y, d;
		x = 0, y = R;
		d = 5 - 4 * R;
		GLint dir[4][2] = {
	   { 1, 1 }, { 1,-1 },
	   { -1,1}, {-1,-1 }
		};		//使用循环来计算对称的区域
		for (GLint i = 0; i < 2; i++) {//这里x=0所以只需要循环2次就行了。
			glVertex2i(Rx, Ry + y * dir[i][1]);
			glVertex2i(Rx + y * dir[i][1], Ry);
		}
		while (x <= y) {
			if (d >= 0) {
				x += 1;
				y -= 1;
				d += 8 * (x - y) + 20;
			}
			else {
				x += 1;
				d += 8 * x + 12;
			}
			for (GLint i = 0; i < 4; i++) {
				glVertex2i(Rx + x * dir[i][0], Ry + y * dir[i][1]);
				glVertex2i(Rx + y * dir[i][1], Ry + x * dir[i][0]);
			}
		}
		glFlush();
		glEnd();
	}
	/*
	PRESENHAM画圆算法
	*/
	/// <summary>
	/// BRESENHAM画圆算法
	/// </summary>
	/// <param name="Rx">圆心</param>
	/// <param name="Ry"></param>
	/// <param name="R ">半径</param>
	void BRESENHAM_Circle(GLint Rx, GLint Ry, GLint R) {
		glBegin(GL_POINTS);
		GLint x, y;
		x = 0, y = R;
		GLint d = 3 - 2 * R;
		GLint dir[4][2] = {
	   { 1, 1 }, { 1,-1 },
	   { -1,1}, {-1,-1 }
		};		//使用循环来计算对称的区域
		for (GLint i = 0; i < 2; i++) {
			glVertex2i(Rx, Ry + y * dir[i][1]);
			glVertex2i(Rx + y * dir[i][0], Ry);
		}
		while (x < y) {
			if (d >= 0) {
				x += 1;
				y -= 1;
				d += 4 * (x - y) + 10;
			}
			else {
				x += 1;
				d += 4 * x + 6;
			}
			for (GLint i = 0; i < 4; i++) {
				glVertex2i(Rx + x * dir[i][0], Ry + y * dir[i][1]);
				glVertex2i(Rx + y * dir[i][0], Ry + x * dir[i][1]);
			}
		}
		glFlush();
		glEnd();
	}
};

//圆弧
class Arc {
public:
	Arc(GLint P1x, GLint P1y, GLint P2x, GLint P2y, GLint P3x, GLint P3y):Rx(P1x),Ry(P1y),start(P2x,P2y){
		R = static_cast<GLint>(sqrt((P1x - P2x) * (P1x - P2x) + (P1y - P2y) * (P1y - P2y)));
		double distance = sqrt((P3x - P1x) * (P3x - P1x) + (P3y - P1y) * (P3y - P1y));
		end.x = static_cast<GLint>(((P3x - P1x) * R) / distance + P1x);
		end.y = static_cast<GLint>(((P3y - P1y) * R) / distance + P1y);
	}
	Arc(GLint rx=0, GLint ry=0, GLint r=0, GLint startx=0, GLint starty=0, GLint endx=0, GLint endy=0)
	:Rx(rx),Ry(ry),R(r),start(startx,starty),end(endx,endy){	}

	void setData(GLint P1x, GLint P1y, GLint P2x, GLint P2y, GLint P3x, GLint P3y) {
		Rx = P1x, Ry = P1y;
		start.x = P2x, start.y = P2y;
		R = static_cast<GLint>(sqrt((P1x - P2x) * (P1x - P2x) + (P1y - P2y) * (P1y - P2y)));
		double distance = sqrt((P3x - P1x) * (P3x - P1x) + (P3y - P1y) * (P3y - P1y));
		end.x = static_cast<GLint>(((P3x - P1x) * R) / distance + P1x);
		end.y = static_cast<GLint>(((P3y - P1y) * R) / distance + P1y);
	}
	void setData(GLint rx, GLint ry, GLint r, GLint startx, GLint starty, GLint endx, GLint endy) {
		Rx = rx, Ry = ry, R = r;
		start.x = startx, start.y = starty;
		end.x = endx, end.y = endy;
	}
	void Draw() {
		DrawCircleArc(Rx, Ry, R, start.x, start.y, end.x, end.y);
	}

private:
	GLint Rx, Ry, R;
	myPoint start, end;

	/*
	圆弧段扫描转换
	*/

	void DrawCircleArc(GLint Rx, GLint Ry, GLint R, GLint startx, GLint starty, GLint endx, GLint endy) {
		//圆心移动到原点，并把y轴反转
		GLint sx = startx - Rx, sy = Ry - starty, ex = endx - Rx, ey = Ry - endy;
		//第一步查找两个端点所在区域编号
		GLint T[2][2], T1[2][2];
		GLint pt0x = 0, pt0y = 0, pt1x = 0, pt1y = 0;
		//找到对应的区间编号（使用移动后的坐标）
		GLint ins = FIndIndexOfArc(sx, sy);
		GLint ine = FIndIndexOfArc(ex, ey);
		//保证结束区域的值大于开始区域，最后使用index%8得到计算区域
		if (ine < ins)ine += 8;
		if (ins == ine) {//同一个区域
			//首先镜像到0号区域，然后镜像到对应的区域
			GetMatrixT(T, ins);				//获得0区域的镜像矩阵
			GetTtoPt(T, sx, sy, pt0x, pt0y, GL_TRUE);
			GetTtoPt(T, ex, ey, pt1x, pt1y, GL_TRUE);
			if (pt0x < pt1x) {//保证了startx <= endx
				SetArc(Rx, Ry, R, 3, T, pt0x, pt0y, pt1x, pt1y);			//画弧线上的一段，从pt0到pt1
			}
			else {
				SetArc(Rx, Ry, R, 3, T, pt1x, pt1y, pt0x, pt0y);
			}
		}
		else {
			//1.起始段画圆弧
			GetMatrixT(T, ins);
			GetTtoPt(T, sx, sy, pt0x, pt0y, GL_TRUE);
			if (!(ins % 2)) {
				SetArc(Rx, Ry, R, 2, T, pt0x, pt0y);			//画弧线上的一段，从(0 ,R )到pt1
			}
			else {
				SetArc(Rx, Ry, R, 1, T, pt0x, pt0y);
			}
			//2.中间画弧段
			for (int i = ins + 1; i < ine; i++) {
				GetMatrixT(T1, i);
				SetArc(Rx, Ry, R, 0, T1);			//画整段弧线
			}
			//3.终止点画弧段
			GetMatrixT(T, ine);
			GetTtoPt(T, ex, ey, pt1x, pt1y, GL_TRUE);
			if (!(ine % 2)) {
				SetArc(Rx, Ry, R, 1, T, pt1x, pt1y);			//画弧线上的一段，从(0 ,R )到pt1
			}
			else {
				SetArc(Rx, Ry, R, 2, T, pt1x, pt1y);
			}
		}
	}

	/// <summary>
	/// 获取对应的区间
	/// </summary>
	GLint FIndIndexOfArc(GLint x, GLint y) {
		if (x > 0 && y >= 0) {
			if (x < y)return 0;
			else return 1;
		}
		else if (x >= 0 && y < 0) {
			if (x > -y)return 2;
			else return 3;
		}
		else if (x < 0 && y <= 0) {
			if (x > y)return 4;
			else return 5;
		}
		else {
			if (-x > y)return 6;
			else return 7;
		}
	}
	/// <summary>
	/// 获取对应区间的变换矩阵
	/// </summary>
	/// <param name="T">得到的变换矩阵</param>
	/// <param name="index">区间id</param>
	void GetMatrixT(GLint T[2][2], GLint index) {
		int i = index % 8;
		if (i == 0) {
			T[0][0] = 1, T[0][1] = 0,
				T[1][0] = 0, T[1][1] = 1;
		}
		else if (i == 1) {
			T[0][0] = 0, T[0][1] = 1,
				T[1][0] = 1, T[1][1] = 0;
		}
		else if (i == 2) {
			T[0][0] = 0, T[0][1] = 1,
				T[1][0] = -1, T[1][1] = 0;
		}
		else if (i == 3) {
			T[0][0] = 1, T[0][1] = 0,
				T[1][0] = 0, T[1][1] = -1;
		}
		else if (i == 4) {
			T[0][0] = -1, T[0][1] = 0,
				T[1][0] = 0, T[1][1] = -1;
		}
		else if (i == 5) {
			T[0][0] = 0, T[0][1] = -1,
				T[1][0] = -1, T[1][1] = 0;
		}
		else if (i == 6) {
			T[0][0] = 0, T[0][1] = -1,
				T[1][0] = 1, T[1][1] = 0;
		}
		else if (i == 7) {
			T[0][0] = -1, T[0][1] = 0,
				T[1][0] = 0, T[1][1] = 1;
		}
	}
	/// <summary>
	/// 进行矩阵乘法
	/// </summary>
	/// <param name="T">镜像矩阵</param>
	/// <param name="sx">开始的位置</param>
	/// <param name="sy"></param>
	/// <param name="ptx">得到的位置</param>
	/// <param name="pty"></param>
	/// <param name="mode">是否进行逆运算</param>
	void GetTtoPt(GLint T[2][2], GLint sx, GLint sy, GLint& retx, GLint& rety, bool mode=GL_FALSE) {
		if (mode) {		//求逆
			//主对角线互换除以行列式，副对角线互换乘-1，除以行列式
			GLint T1[2][2], temp;
			T1[0][0] = T[1][1];
			T1[1][1] = T[0][0];
			T1[1][0] = T[1][0] * -1;
			T1[0][1] = T[0][1] * -1;
			temp = T[0][0] * T[1][1] - T[1][0] * T[0][1];
			retx = (sx * T1[0][0] + sy * T1[0][1]) / temp;
			rety = (sx * T1[1][0] + sy * T1[1][1]) / temp;
		}
		else {
			retx = (sx * T[0][0] + sy * T[0][1]);
			rety = (sx * T[1][0] + sy * T[1][1]);
		}
	}
	/// <summary>
	/// 绘制弧段
	///0 画1/8圆  是连续段的中部
	///1 画0到start圆 是奇数段的开始   
	///2 画start到1/8x圆  是偶数段的开始
	///3 画start到end圆   是连续段的全部
	/// </summary>
	void SetArc(GLint Rx, GLint Ry, GLint R, GLint ArcType, GLint T[2][2], GLint startx=0, GLint starty=0, GLint endx=0, GLint endy=0) {
		GLint x, y, d;
		GLint p0x, p0y, p1x, p1y;
		glBegin(GL_POINTS);

		//0 画1/8圆  是连续段的中部
		//1 画0到start圆 是奇数段的开始   
		//2 画start到1/8x圆  是偶数段的开始
		//3 画start到end圆   是连续段的全部
		if (ArcType == 0 || ArcType == 1) {
			x = 0;
			y = R;
			d = 5 - 4 * R;
			glColor3f(1.0, 0, 0);
		}
		else {
			glColor3f(0, 1.0, 0);
			x = startx;
			y = starty;
			d = 8 * x - 4 * y + 5;
		}
		p0x = x, p0y = y;
		GetTtoPt(T, p0x, p0y, p1x, p1y);
		glVertex2i(Rx + p1x, Ry - p1y);
		while (1) {
			if (ArcType == 0 || ArcType == 2) {		//生成到镜像线八分之一处。
				if (x > y)break;
			}
			else if (ArcType == 1) {		//生成（0,R）到start的圆弧
				glColor3f(0, 1.0, 0);
				if (x > startx)break;
			}
			else if (ArcType == 3) {						//生成到end的圆弧
				glColor3f(0, 0, 1.0);
				if (x > endx)break;
			}
			if (d >= 0) {
				x += 1, y -= 1;
				d += 8 * (x - y) + 20;
			}
			else {
				x += 1;
				d += 8 * x + 12;
			}
			p0x = x;
			p0y = y;
			GetTtoPt(T, p0x, p0y, p1x, p1y);
			glVertex2i(Rx + p1x, Ry - p1y);
		}
		glColor3f(1.0, 0, 0);
		glFlush();
		glEnd();
	}
};

//椭圆
class Ellipt {
public:
	Ellipt(GLint P1x=0,GLint P1y=0,GLint P2x=0,GLint P2y=0)
	{
		m_a = abs(static_cast<GLint>(P1x - P2x)/2);
		m_b= abs(static_cast<GLint>(P1y - P2y)/2);
		//四舍五入
		m_R.x = static_cast<GLint>(((P1x + P2x) / 2.0 + 0.5));
		m_R.y = static_cast<GLint>(((P1y + P2y) / 2.0 + 0.5));
	}
	void setData(GLint P1x, GLint P1y, GLint P2x, GLint P2y) {
		m_a = abs(static_cast<GLint>(P1x - P2x)/2);
		m_b = abs(static_cast<GLint>(P1y - P2y)/2);
		//四舍五入
		m_R.x = static_cast<GLint>(((P1x + P2x) / 2.0 + 0.5));
		m_R.y = static_cast<GLint>(((P1y + P2y) / 2.0 + 0.5));
	}
	void Draw(){
		MidPt_Elliptse(m_R, m_a, m_b);
	}
private:
	myPoint m_R;
	GLint m_a,m_b;
	void MidPt_Elliptse(myPoint cPt, GLint a, GLint b) {
		glBegin(GL_POINTS);
		GLint x, y, temp = a * 7 / 10, dir[4][2] = {
		{ 1, 1 }, { 1,-1 },
		{ -1,1}, {-1,-1 }
		};
		double d;
		x = 0, y = b;
		d = b * b + a * a * (-b + 0.25);
		for (int i = 0; i < 4; i++) {
			glVertex2i(cPt.x + x * dir[i][0], cPt.y + y * dir[i][1]);
		}
		while ((b * b * (x + 1)<a * a * (y + 0.5))) {
			if (d > 0)
			{
				x += 1;
				y -= 1;
				d += b * b * (2 * x + 3) + a * a * (-2 * y + 2);
			}
			else {
				x += 1;
				d += b * b * (2 * x + 3);
			}
			for (int i = 0; i < 4; i++) {
				glVertex2i(cPt.x + x * dir[i][0], cPt.y + y * dir[i][1]);
			}
		}
		while (y > 0) {
			if (d >= 0) {
				y -= 1;
				d += a * a * (-2 * y + 3);
			}
			else {
				x += 1;
				y -= 1;
				d += a * a * (-2 * y + 3) + b * b * (2 * x + 2);
			}
			for (int i = 0; i < 4; i++) {
				glVertex2i(cPt.x + x * dir[i][0], cPt.y + y * dir[i][1]);
			}
		}
		glFlush();
		glEnd();
	}
};
#endif