#ifndef POINT_HPP
#define POINT_HPP
#include <GL/freeglut.h>

class myPoint;

//A point or a vector in R^2
class myPoint {
public:
	myPoint(GLint X=0,GLint Y=0 ):x(X),y(Y){}
	GLint x, y;
};
myPoint operator+(const myPoint& a, const myPoint& b) {
	return myPoint(a.x + b.x, a.y + b.y);
}
myPoint operator-(const myPoint& a, const myPoint& b) {
	return myPoint(a.x - b.x, a.y - b.y);
}
GLint operator*(const myPoint& a, const myPoint& b) {
	return a.x * b.x + a.y * b.y;
}
#endif