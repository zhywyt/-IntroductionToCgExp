#pragma once
#include "Cell.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <time.h>
#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

namespace MeshKernel {
	// 网格的概念
	class Mesh {
	public:
		iGameVertex BBoxMin, BBoxMax;
		iGameEdgeHandle getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2) {// return the edge whose vex is vh1 and vh2
			iGameEdgeHandle ret(-1);
			for (auto& eh : NeighborEh(vh1)) {
				auto& e = edges_[eh];
				auto adj = e.vh1() + e.vh2() - vh1;
				if (vh2 == adj) {
					ret = eh;
					break;
				}
			}
			return ret;
		}

		void initBBox();
		inline bool isValid(iGameVertexHandle _vh) { return vertices_.count(_vh); }
		inline bool isValid(iGameEdgeHandle _eh) { return edges_.count(_eh); }
		inline bool isValid(iGameFaceHandle _fh) { return faces_.count(_fh); }
		inline bool isSingular(iGameVertexHandle vh) {
			return NeighborEhOfVertex_[vh].empty();
		}


		/*=========================读写元素===============================*/
			// 读取ID为i的顶点
		iGameVertex& vertices(iGameVertexHandle _vh);
		const iGameVertex vertices(iGameVertexHandle _vh) const;        // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用
		// 读取ID为i的边
		iGameEdge& edges(iGameEdgeHandle _eh);
		const iGameEdge& edges(iGameEdgeHandle _eh) const;
        bool edges_vaild(iGameEdgeHandle _eh);
		// 读取ID为i的面
		iGameFace& faces(iGameFaceHandle _fh);
		const iGameFace faces(iGameFaceHandle _fh) const;

        bool faces_vaild(iGameFaceHandle _fh);

		size_t vsize() const { return vertices_.size(); }
		size_t esize() const { return edges_.size(); }
		size_t fsize() const { return faces_.size(); }
		const std::unordered_map<iGameVertexHandle, iGameVertex>& allvertices() const { return vertices_; }
		const std::unordered_map<iGameEdgeHandle, iGameEdge>& alledges() const { return edges_; }
		const std::unordered_map<iGameFaceHandle, iGameFace>& allfaces() const { return faces_; }

		/*====================根据元素得到对应ID=========================*/
		const iGameVertexHandle vertexhandle(iGameVertex _vertex) const;
		const iGameEdgeHandle edgehandle(iGameEdge& _edge) const;
		const iGameFaceHandle facehandle(iGameFace& _face) const;
		/*======================得到邻接关系============================*/
		// 顶点的邻接点
		std::unordered_set<iGameVertexHandle> NeighborVh(iGameVertexHandle _vh);
		// 顶点的邻接边
		std::unordered_set<iGameEdgeHandle>& NeighborEh(iGameVertexHandle _vh);
		// 顶点的邻接面
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameVertexHandle _vh);
		// 根据顶点的一条边得到对边的点： ·—·
		iGameVertexHandle NeighborVhFromEdge(iGameVertexHandle _vh, iGameEdgeHandle _eh);
        // 根据顶点的一条边得到对点的边： —·— ，四边形专用
        iGameEdgeHandle NeighborEhFromVertex(iGameEdgeHandle _eh, iGameVertexHandle _vh);
		// 边的邻接边
		std::unordered_set<iGameEdgeHandle> NeighborEh(iGameEdgeHandle _eh);
		// 边的邻接面
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameEdgeHandle _eh);
		// 面的邻接面
		std::unordered_set<iGameFaceHandle> NeighborFh(iGameFaceHandle _fh);// share common edge
		std::unordered_set<iGameFaceHandle> Neighbor2Fh(iGameFaceHandle _fh);// share common vertex
        // 根据面的一条边得到对边的面： 口|口
        iGameFaceHandle NeighborFhFromEdge(iGameFaceHandle _fh, iGameEdgeHandle _eh);
        // 根据面的一条边得到对边： |口| ，四边形专用，三角形或多边形会返回随机的另一条边
        iGameEdgeHandle OppositeEhFromEdge(iGameFaceHandle _fh, iGameEdgeHandle _eh);
        // 根据两个顶点找到边，没找到就返回-1
        iGameEdgeHandle getEhFromTwoVh(iGameVertexHandle _vh1, iGameVertexHandle _vh2);

        // 获得边的向量（未归一化）
        iGameVertex getEdgeVector(iGameEdge _e);

		/*=========================添加元素=============================*/
		iGameVertexHandle AddVertex(const iGameVertex& _v);
		iGameEdgeHandle AddEdge(const iGameVertexHandle& _vh1, const iGameVertexHandle& _vh2);
		iGameFaceHandle AddFace(const std::vector<iGameVertexHandle>& _vhs);

		/*=========================删除元素=============================*/
		// 删除低级元素时删除一切邻接的高级元素
		// 删除高级元素时保留仍被其他元素使用的低级元素
		iGameVertexHandle DeleteVertex(iGameVertexHandle _vh);
		iGameEdgeHandle DeleteEdge(iGameEdgeHandle _eh);
		iGameFaceHandle DeleteFace(iGameFaceHandle _fh);

		/*=========================生成唯一ID========================*/
		iGameVertexHandle GenVertexHandle() { return (iGameVertexHandle)VertexHandleID_++; }
		iGameEdgeHandle GenEdgeHandle() { return (iGameEdgeHandle)EdgeHandleID_++; }
		iGameFaceHandle GenFaceHandle() { return (iGameFaceHandle)FaceHandleID_++; }

		size_t VertexSize() { return vertices_.size(); }
		size_t EdgeSize() { return edges_.size(); }
		size_t FaceSize() { return faces_.size(); }
	protected:
		/*============下一个可使用的ID=========*/
		int VertexHandleID_ = 0;
		int EdgeHandleID_ = 0;
		int FaceHandleID_ = 0;
		std::unordered_set<iGameVertexHandle> empty_vhs;
		std::unordered_set<iGameEdgeHandle> empty_ehs;
		std::unordered_set<iGameFaceHandle> empty_fhs;

		/*=============修改邻接关系============*/
		// 将该面添加至面所包含的点和边的相邻面中
		void AddFace2Neighbor(const iGameFaceHandle& _fh);
		// 将该边添加至边所包含的点相邻边中
		void AddEdge2Neighbor(const iGameEdgeHandle& _eh);
		// Todo: Delete Neighbor
		void DeleteFace2Neighbor(const iGameFaceHandle& _fh);
		void DeleteEdge2Neighbor(const iGameEdgeHandle& _eh);
	protected:
		// handle到元素的对应
		std::unordered_map<iGameVertexHandle, iGameVertex> vertices_;
		std::unordered_map<iGameEdgeHandle, iGameEdge> edges_;
		std::unordered_map<iGameFaceHandle, iGameFace> faces_;

		// 元素到handle的对应
		std::unordered_map<iGameVertex, iGameVertexHandle> Vertex2Vh_;
		std::unordered_map<iGameEdge, iGameEdgeHandle> Edge2Eh_;
		std::unordered_map<iGameFace, iGameFaceHandle> Face2Fh_;

		// 邻接关系
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameEdgeHandle>> NeighborEhOfVertex_;          //点的邻接边
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameFaceHandle>> NeighborFhOfVertex_;          //点的邻接面
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameFaceHandle>> NeighborFhOfEdge_;              //边的邻接面
	protected:
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) = 0;                              // 必须重写
		Mesh& operator=(const Mesh& _mesh);


    };

	// 曲面网格
	class SurfaceMesh : public Mesh {
	public:
		/*==========================构造函数==============================*/
		SurfaceMesh() {};
		SurfaceMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _faces) {
			InitMesh(_vertices, _faces);
		}
		//SurfaceMesh(const SurfaceMesh& _surfacemesh);
		/*=============初始化网格=============*/
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) override;
		SurfaceMesh& operator=(const SurfaceMesh& _surfacemesh);

		bool isOnBoundary(iGameEdgeHandle);
		bool isOnBoundary(iGameVertexHandle);
		bool isOnBoundary(iGameFaceHandle);

		void genNormal(iGameFaceHandle);
		void genNormal(iGameVertexHandle);
		void genAllFacesNormal();
		void genAllVerticesNormal();// 自然会生成所有面的法向量

		bool isClosure();

		void genAllEdgesLength();
		void genLength(iGameEdgeHandle);
		double getLength(iGameEdgeHandle);
        //仅支持计算三角形和四边形的面积
        double getFaceArea(iGameFaceHandle fh);
		void addNoise(double);

		iGameEdgeHandle getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2);

        std::unordered_map<int, int> updateAllHandles();

		bool isConnected(iGameVertexHandle, iGameVertexHandle);
		bool isConnected(iGameEdgeHandle, iGameEdgeHandle);

		bool isTriangleMesh();
		size_t getBoundaryVerticesCount();
	};

	// 体网格
	class VolumeMesh : public Mesh {
	public:
		/*==========================构造函数==============================*/
		VolumeMesh() {};
		VolumeMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _cells) {
			InitMesh(_vertices, _cells);
			//std::cout << "InitMesh(_vertices, _cells) : debug " <<std::endl;
			for(auto i : allcells()){
			    for(auto j : i.second.getVertexHandle()){
                    assert(j.idx() < VertexSize());
			    }
			}
		}
		void updateAllHandles();
        void updateAllHandles_2(std::vector<int> &vh_list);
		/*=============初始化网格=============*/
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) override;
		VolumeMesh& operator=(const VolumeMesh& _volumemesh);

	public:
		iGameCell& cells(iGameCellHandle _ch);
		const iGameCell cells(iGameCellHandle _ch) const;
		inline bool isValid(iGameCellHandle _ch) { return cells_.count(_ch); }

		size_t CellSize() const { return cells_.size(); }
		const std::unordered_map <iGameCellHandle, iGameCell>& allcells() const { return cells_; }

		const iGameCellHandle cellhandle(iGameCell& _cell) const;

		std::unordered_set<iGameCellHandle> NeighborCh(iGameVertexHandle _vh);// 点邻接体
		std::unordered_set<iGameCellHandle> NeighborCh(iGameEdgeHandle _eh);// 边邻接体
		std::unordered_set<iGameCellHandle> NeighborCh(iGameFaceHandle _fh);// 面邻接体
		std::unordered_set<iGameCellHandle> NeighborCh(iGameCellHandle _ch);// 体邻接体

		iGameCellHandle AddCell(const std::vector<iGameVertexHandle>& _vhs);// 加体
		iGameCellHandle AddCell(const std::vector< std::vector<iGameVertexHandle>>& _vhs);

		iGameVertexHandle DeleteVertex(const iGameVertexHandle& _vh);// 删点
		iGameEdgeHandle DeleteEdge(const iGameEdgeHandle& _eh);// 删边
		iGameFaceHandle DeleteFace(const iGameFaceHandle& _fh);// 删面
		iGameCellHandle DeleteCell(const iGameCellHandle& _ch);// 删体

		iGameCellHandle GenCellHandle() { return (iGameCellHandle)CellHandleID_++; }

		bool isConnected(iGameFaceHandle& fh1, iGameFaceHandle& fh2);
		bool isConnected(iGameEdgeHandle& eh1, iGameEdgeHandle& eh2);
		bool isConnected(iGameVertexHandle& vh1, iGameVertexHandle& vh2);
		bool isConnected(const iGameVertexHandle& vh1, const iGameVertexHandle& vh2);

		bool isOnBoundary(iGameCellHandle ch);
		bool isOnBoundary(iGameFaceHandle fh);
		bool isOnBoundary(iGameEdgeHandle eh);
		bool isOnBoundary(iGameVertexHandle vh);

		iGameVertex getCenter(iGameCellHandle ch);
		iGameVertex getQuadNormal(iGameFaceHandle fh);
		double getQuadArea(iGameFaceHandle fh);
		double getLength(iGameEdgeHandle eh);

        void genNormal(iGameFaceHandle);
		void genNormal(iGameVertexHandle);
		void genAllFacesNormal();
		void genAllVerticesNormal();

        size_t CellSize() { return cells_.size(); }

	protected:
		int CellHandleID_ = 0;
		std::unordered_set<iGameCellHandle> empty_chs;
		std::vector<std::vector<iGameVertexHandle>> surface_faces_;
		void AddCell2Neighbor(const iGameCellHandle& _ch);
		void DeleteCell2Neighbor(const iGameCellHandle& _ch);

		std::unordered_map<iGameCellHandle, iGameCell> cells_;
		std::unordered_map<iGameCell, iGameCellHandle> Cell2Ch_;

		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameCellHandle>> NeighborChOfVertex_;          //点的邻接体
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameCellHandle>> NeighborChOfEdge_;              //边的邻接体
		std::unordered_map<iGameFaceHandle, std::unordered_set<iGameCellHandle>> NeighborChOfFace_;              //面的邻接体
        iGameFaceHandle check2();
    };

	class TriMesh : SurfaceMesh {};
	class QuadMesh : SurfaceMesh {};

    class TetMesh : VolumeMesh {};
//	class VolumeMesh : public VolumeMesh {
//	public:
//        VolumeMesh() {};
//		VolumeMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _cells) {
//			InitMesh(_vertices, _cells);
//			/*InitHedra(_vertices, _cells, _surface_faces);*/
//			printf("init tetrahedron success\n");
//		}
//
//		/*void InitHedra(const std::vector<iGameVertex>& _vertices, std::vector<std::vector<std::vector<iGameVertexHandle>>>& _cells,
//			const std::vector<std::vector<iGameVertexHandle>>& _surface_faces);*/
//		void genNormal(iGameFaceHandle);
//		void genNormal(iGameVertexHandle);
//		void genAllFacesNormal();
//		void genAllVerticesNormal();
//
//	};


	class HexMesh : VolumeMesh {};
}
