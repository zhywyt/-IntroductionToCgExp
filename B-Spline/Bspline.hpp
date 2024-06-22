#ifndef BSPLINE_HPP
#define BSPLINE_HPP
#include <vector>
#include "Point.hpp"
#include "DrawLine.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

class BSpline{
    public:
        friend ostream& operator<<(ostream&,BSpline&);

        
        BSpline( vector<myPoint> points=vector<myPoint>(),int n=-1, int p=3, int head=2, int node=1000)
        :m_p(p),m_head(head),m_node(node),m_points(points),m_n(n)
        {
            m_nodeVec = getVec(m_n,m_p,m_head);
            m_n=points.size()-1;
        }
        myPoint& operator [](int index){
            return m_points[index];
        }
        GLint getPoint(GLint x,GLint y){
            for(int i=0;i<m_points.size();i++){
                cout<<"check x:"<<m_points[i].x<<" check y:"<<m_points[i].y<<endl;
                if((m_points[i].x-x)*(m_points[i].x-x)+(m_points[i].y-y)*(m_points[i].y-y)<200)return i;
            }
            return -1;
        }
        bool save(string filename){
            ofstream os(filename);
            if(!os)return false;
            else{
                os<<*this;
            }
            /*
            int m_n, m_p, m_head, m_node;
            vector<GLdouble> m_nodeVec;
            vector<myPoint> m_points;
            */
            return true;
        }
        bool load(string filename){
            ifstream is(filename);
            if(!is)return false;
            vector<myPoint>points;
            GLint n,p,head,node,x,y;
            is>>n>>p>>head>>node;
            while(is>>x>>y){
                points.push_back(myPoint(x,y));
            }
            if(points.size()-1!=n)return false;
            BSpline temp(points,n,p,head,node);
            swap(*this,temp);
            return true;
        }
        void erase(GLint index){
            if(index>=0&&index<m_points.size()){
                m_points.erase(m_points.begin()+index);
                m_n--;
                m_node=m_n*200;
            }
        }
        void pushControlP(myPoint p){
            m_points.push_back(p);
            m_n++;
            m_node=m_n*200;
        }        
        void clear(){
            m_points.clear();
            m_nodeVec.clear();
            m_n=-1;
            m_head=2;
            m_p=3;
            m_node=1000;
        }
        void draw(int p=-1, int head=-1, int node = -1)
        {
            if(p==-1)p=m_p;
            if(head==-1)head=m_head;
            if(node==-1)node=m_node;
            glPointSize(10);
            glBegin(GL_POINTS);
            for(int i=0;i<m_points.size();i++){
                glVertex2i(m_points[i].x,m_points[i].y);
            }
            glEnd();
            glColor3f(0.0f,1.0f,0.0f);
            glPointSize(4);
            for(int i=0,j=1;j<m_points.size();i++,j++){
                BRESENHAM_Line(m_points[i].x,m_points[i].y,m_points[j].x,m_points[j].y);
            }
            glColor3f(0.0f,0.0f,0.0f);
            if(head*2 > m_n + m_p +1)return;       //点不够
            if(p > m_n)return;          //点不够


            m_nodeVec = getVec(m_n, m_p, head);
            
            GLdouble basis;
            vector<GLdouble> rx(node, 0), ry(node, 0);
            // GLdouble u_begin = nodeVec[p - 1];
            // GLdouble u_end = nodeVec[n + 1];

            GLdouble u_begin = m_nodeVec[p-1];
            GLdouble u_end = m_nodeVec[m_n+1];
            GLdouble dis = (u_end - u_begin) / node;
            GLdouble u = u_begin;
            for (int i = 0; i < node; i++)
            {
                for (int j = 0; j <= m_n; j++)
                {
                    basis = b_spline_basis(j, p, u, m_nodeVec);
                    // cout << "i: " << i << " j: " << j << " basis: " << basis << endl;
                    rx[i] += m_points[j].x * basis;
                    ry[i] += m_points[j].y * basis;
                }
                u += dis;
                // cout << "x[u]: " << rx[i] << " y[u]: " << ry[i] << endl;
            }
            glBegin(GL_POINTS);
            for (int i = 0; i < node; i++){
                glVertex2d(rx[i], ry[i]);
            }
            glEnd();
            glFlush();
        }
        bool setPara(char c,int data){
            switch (c)
            {
            case 'p':
            case 'P':
                m_p = data;
                m_nodeVec =  getVec(m_n,m_p,m_head);
                return true;
                
            case 'h':
            case 'H':
                m_head = data;
                m_nodeVec =  getVec(m_n,m_p,m_head);
                return true;
            default:
                break;
            }
            return false;
        }
        void showInfo(ostream&os){
            os<<"n\t: "<<m_n+1<<"\t\t\tnumber of control points.\n"
                <<"p\t: "<<m_p<<"\t\t\trank of BSpline.\n"
                <<"h\t: "<<m_head<<"\t\t\trepeat of head and tail.\n"
                <<"nv\t: "<<m_nodeVec.size()<<"\t\t\tnumber of node vector.\n"
                <<"head vector:\n";
            for(auto i :m_nodeVec){
                os<<i<<" ";
            }
            os<<endl;

        }



    // private:
        static  myPoint FALSEPOINT;
        vector<GLdouble> getVec(int n, int p, int head){
            if(head*2 > m_n + m_p + 1)return vector<GLdouble>();        //点不够
            if(p > m_n)return vector<GLdouble>();                       //点不够

            vector<GLdouble> ans;
            GLdouble u_temp = 0;
            GLdouble dis = 1.0 / (n + p - 2*(head-1) + 1 );
            //head 
            for (int i = 0; i < head; i++)
            {
                ans.push_back(0);
            }
            //n + p + 1 - 2*(head )
            u_temp=dis;
            for (int i = 0; i < n + p + 1 - 2 * (head); i++)
            {
                ans.push_back(u_temp);
                u_temp += dis;
            }
            // cout<<"u_temp: "<<u_temp<<endl;
            for (int i = 0; i < head; i++)
            {
                ans.push_back(1);
            }
            //n + p + 1
            return ans;
        }

        GLdouble b_spline_basis(int i, int p, GLdouble u, const vector<GLdouble> &nodeVec){
            // cout<<"u: "<<u<<endl;
            if (p == 1){
                if (nodeVec[i] <= u && u <= nodeVec[i + 1])
                    return 1;
                else
                    return 0;
            }
            if (u < nodeVec[i] || u > nodeVec[i + p])
                return 0;
            GLdouble len1 = nodeVec[i + p - 1] - nodeVec[i];
            GLdouble len2 = nodeVec[i + p] - nodeVec[i + 1];
            GLdouble alpha, beta;
            if (len1 <= 1e-10){
                alpha = 0;
            }
            else{
                alpha = (u - nodeVec[i]) / len1;
            }
            if (len2 <= 1e-10){
                beta = 0;
            }
            else{
                beta = (nodeVec[i + p] - u) / len2;
            }
            return (alpha)*b_spline_basis(i, p - 1, u, nodeVec) + (beta)*b_spline_basis(i + 1, p - 1, u, nodeVec);
        }
         /*
        点的数量-1
        曲线的阶数
        头节点重叠数
        绘制曲线的点数
        */
        int m_n, m_p, m_head, m_node;
        vector<GLdouble> m_nodeVec;
        vector<myPoint> m_points;
};

myPoint BSpline::FALSEPOINT{-1,-1};
ostream& operator<<(ostream&os,BSpline&bspline){
    if(!os)return os;
    bspline.m_n=bspline.m_points.size()-1;
    os<<bspline.m_n<<" "<<bspline.m_p<<" "<<bspline.m_head<<" "<<bspline.m_node<<endl;
    for(auto& i :bspline.m_points){
        os<<i.x<<" "<<i.y<<"\n";
    }
    os<<endl;
    return os;
}


#endif
