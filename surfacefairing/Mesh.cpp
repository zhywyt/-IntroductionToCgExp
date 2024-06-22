#include"Mesh.h"
#include"string"
#include <queue>

// Mesh 定义
namespace MeshKernel {

    void Mesh::initBBox() {
        double bbox_min_x = 99999999, bbox_min_y = 99999999, bbox_min_z = 99999999;
        double bbox_max_x = -99999999, bbox_max_y = -99999999, bbox_max_z = -99999999;
        for (auto& vp : vertices_) {
            bbox_min_x = std::min(bbox_min_x, vp.second.x());
            bbox_min_y = std::min(bbox_min_y, vp.second.y());
            bbox_min_z = std::min(bbox_min_z, vp.second.z());
            bbox_max_x = std::max(bbox_max_x, vp.second.x());
            bbox_max_y = std::max(bbox_max_y, vp.second.y());
            bbox_max_z = std::max(bbox_max_z, vp.second.z());
        }
        BBoxMin = iGameVertex(bbox_min_x, bbox_min_y, bbox_min_z);
        BBoxMax = iGameVertex(bbox_max_x, bbox_max_y, bbox_max_z);
        printf("volume mesh: BBox: (%.3f, %.3f, %.3f) --> (%.3f, %.3f, %.3f)\n",
               BBoxMin.x(), BBoxMin.y(), BBoxMin.z(), BBoxMax.x(), BBoxMax.y(), BBoxMax.z());
    }

    iGameFaceHandle Mesh::AddFace(const std::vector<iGameVertexHandle>& _vhs) {
        //std::cout << " 现在加的面中的点为 : " << std::endl;
        //for (auto& vh : _vhs) {
        //	std::cout << vh << " : ";
        //}
        //std::cout << std::endl;
        std::vector<iGameEdgeHandle> ehs(_vhs.size());
        // 得到该面边的handle
        for (int i = 0; i < _vhs.size(); ++i) {
            if (i == 0) {
                ehs[i] = AddEdge(_vhs[_vhs.size() - 1], _vhs[i]);
            }
            else {
                ehs[i] = AddEdge(_vhs[i], _vhs[i - 1]);
            }
        }
        iGameFace f(_vhs, ehs);
        // 如果该面已经存在，则返回面的handle
        if (Face2Fh_.count(f)) {
            iGameFaceHandle fh = Face2Fh_[f];
//            std::cout << "该面已经存在 序号为 : " << Face2Fh_[f] << "以及面中的各个边的序号为 : " << std::endl;
//            for (int i = 0; i < ehs.size(); i++) {
//            	std::cout << ehs[i] << " ";
//            }
//            std::cout << " | ";
//            for (int i = 0; i < faces_[fh].size(); i++) {
//            	std::cout << faces_[fh].eh(i) << " ";
//            }
//            std::cout << std::endl;
            // 仍然更新 handle 和 面直接的映射
            faces_[fh] = f;											// 建立handle与该面之间的映射
            return Face2Fh_[f];
        }
            // 否则
        else {
            iGameFaceHandle fh = GenFaceHandle();						// 生成一个新的handle
            //std::cout << " 新的面的序号为 : " << fh << " 和面中的各个边的序号为 : " << std::endl;
            //for (int i = 0; i < ehs.size(); i++) {
            //	std::cout << ehs[i] << " ";
            //}
            faces_[fh] = f;											// 建立handle与该面之间的映射
            Face2Fh_[f] = fh;										// 建立该面与handle之间的映射
            AddFace2Neighbor(fh);									// 将该面添加至面所包含的点和边的相邻面中
            //for (int i = 0; i < 4; i++) {
            //	std::cout << faces_[fh].eh(i) << " ";
            //}
            //std::cout << std::endl;
            return Face2Fh_[f];										// 返回该面的handle
        }
    }
    iGameEdgeHandle Mesh::AddEdge(const iGameVertexHandle& vh1, const iGameVertexHandle& vh2) {
        iGameEdge e(vh1, vh2);
        //std::cout << "边的两个端点的序号为 : " << vh1 << " : " << vh2 << std::endl;
        // 如果该边已经存在，则返回该边的handle
        if (Edge2Eh_.count(e)) {
            //std::cout << "此时该边已经存在 且边的序号为 : " << Edge2Eh_[e] << std::endl;
            return Edge2Eh_[e];
        }
            // 否则
        else {
            iGameEdgeHandle eh = GenEdgeHandle();						// 生成一个新的handle
            //std::cout << "生成边的序号为 : " << eh << std::endl;
            edges_[eh] = e;											// 建立handle与该边之间的映射
            Edge2Eh_[e] = eh;										// 建立该边与handle之间的映射
            AddEdge2Neighbor(eh);									// 将该边记录为其两个顶点的邻接边
            return Edge2Eh_[e];										// 返回该边的handle
        }
    }
    iGameVertexHandle Mesh::AddVertex(const iGameVertex& _v) {
        if (Vertex2Vh_.count(_v)) return Vertex2Vh_[_v];
        else {
            iGameVertexHandle vh = GenVertexHandle();
            vertices_[vh] = _v;
            Vertex2Vh_[_v] = vh;
            return Vertex2Vh_[_v];
        }
    }

    iGameVertexHandle Mesh::DeleteVertex(iGameVertexHandle _vh) {
        if (!vertices_.count(_vh)) return iGameVertexHandle(-1);        // 如果该点不存在，则返回-1
        else {
            // 删除相邻边元素（删除边时自然删除了相邻面）
            auto ve = NeighborEhOfVertex_[_vh];
            for (iGameEdgeHandle eh : ve) {
                DeleteEdge(eh);
            }
            // 删除顶点元素和以其为根据的所有邻接关系
            Vertex2Vh_.erase(vertices_[_vh]);
            vertices_.erase(_vh);
            NeighborEhOfVertex_.erase(_vh);
            NeighborFhOfVertex_.erase(_vh);
            return _vh;
        }

    }
    iGameEdgeHandle Mesh::DeleteEdge(iGameEdgeHandle _eh) {
        if (!edges_.count(_eh)) return iGameEdgeHandle(-1);             // 如果该边不存在，则返回-1
        else {
            // 删除邻接关系
            iGameEdge e(edges_[_eh]);
            for (int i = 0; i < 2; ++i) {
                iGameVertexHandle ev = e.vh(i);
                NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
            }
            // 删除相邻面元素
            auto ef = NeighborFhOfEdge_[_eh];
            for (iGameFaceHandle fh : ef) {
                DeleteFace(fh);
            }
            // 删除边元素和以其为根据的所有邻接关系
            Edge2Eh_.erase(edges_[_eh]);
            edges_.erase(_eh);
            NeighborFhOfEdge_.erase(_eh);
            return _eh;
        }
    }
    iGameFaceHandle Mesh::DeleteFace(iGameFaceHandle _fh) {
        if (!faces_.count(_fh)) return iGameFaceHandle(-1);             // 如果该面不存在，则返回-1
        else {                                                     // 如果该面存在，则返回删除的这个面的handle
            // 删除邻接关系
            iGameFace f(faces_[_fh]);
            for (int i = 0; i < f.size(); ++i) {
                iGameVertexHandle fv = f.vh(i);
                iGameEdgeHandle fe = f.eh(i);
                NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
                NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
            }
            // 删除面元素
            Face2Fh_.erase(faces_[_fh]);
            faces_.erase(_fh);
            return _fh;
        }
    }

    void Mesh::AddFace2Neighbor(const iGameFaceHandle& _fh)
    {
        iGameFace f = faces_[_fh];
        size_t n = f.size();
        for (int i = 0; i < n; ++i) {
            NeighborFhOfVertex_[f.vh(i)].insert(_fh);
        }
        for (int i = 0; i < n; ++i) {
            NeighborFhOfEdge_[f.eh(i)].insert(_fh);
        }
    }

    void Mesh::AddEdge2Neighbor(const iGameEdgeHandle& _eh)
    {
        iGameEdge e = edges_[_eh];
        NeighborEhOfVertex_[e.vh1()].insert(_eh);
        NeighborEhOfVertex_[e.vh2()].insert(_eh);

    }

    void Mesh::DeleteFace2Neighbor(const iGameFaceHandle& _fh) {
        iGameFace f = faces_[_fh];
        size_t n = f.size();
        for (int i = 0; i < n; ++i) {
            NeighborFhOfVertex_[f.vh(i)].erase(_fh);
        }
        for (int i = 0; i < n; ++i) {
            NeighborFhOfEdge_[f.eh(i)].erase(_fh);
        }
    }
    void Mesh::DeleteEdge2Neighbor(const iGameEdgeHandle& _eh) {
        iGameEdge e = edges_[_eh];
        NeighborEhOfVertex_[e.vh1()].erase(_eh);
        NeighborEhOfVertex_[e.vh2()].erase(_eh);
    }

    Mesh& Mesh::operator=(const Mesh& _surfacemesh)
    {
        vertices_ = _surfacemesh.vertices_;
        edges_ = _surfacemesh.edges_;
        faces_ = _surfacemesh.faces_;
        Vertex2Vh_ = _surfacemesh.Vertex2Vh_;
        Edge2Eh_ = _surfacemesh.Edge2Eh_;
        Face2Fh_ = _surfacemesh.Face2Fh_;
        NeighborEhOfVertex_ = _surfacemesh.NeighborEhOfVertex_;
        NeighborFhOfVertex_ = _surfacemesh.NeighborFhOfVertex_;
        NeighborFhOfEdge_ = _surfacemesh.NeighborFhOfEdge_;
        VertexHandleID_ = _surfacemesh.VertexHandleID_;
        EdgeHandleID_ = _surfacemesh.EdgeHandleID_;
        FaceHandleID_ = _surfacemesh.FaceHandleID_;
        return *this;
    }

    /*=========================读写元素===============================*/
    // 读取ID为i的顶点
    iGameVertex& Mesh::vertices(iGameVertexHandle _vh) {
        assert(vertices_.count(_vh));
        return vertices_[_vh];
    }
    const iGameVertex Mesh::vertices(iGameVertexHandle _vh) const {
        assert(vertices_.count(_vh));
        return vertices_.find(_vh)->second;                // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用
    }

    // 读取ID为i的边
    iGameEdge& Mesh::edges(iGameEdgeHandle _eh) {
        assert(edges_.count(_eh));
        return edges_[_eh];
    }
    bool Mesh::edges_vaild(iGameEdgeHandle _eh) {
        return edges_.count(_eh);
    }
    const iGameEdge& Mesh::edges(iGameEdgeHandle _eh) const {
        assert(edges_.count(_eh));
        return edges_.find(_eh)->second;
    }
    // 读取ID为i的面
    iGameFace& Mesh::faces(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));
        return faces_[_fh];
    }
    const iGameFace Mesh::faces(iGameFaceHandle _fh) const {
        assert(faces_.count(_fh));
        return faces_.find(_fh)->second;
    }
    bool Mesh::faces_vaild(iGameFaceHandle _fh) {
        return  faces_.count(_fh);
    }
    /*====================根据元素得到对应ID=========================*/
    const iGameVertexHandle Mesh::vertexhandle(iGameVertex _vertex) const {
        if (Vertex2Vh_.find(_vertex) != Vertex2Vh_.end()) return Vertex2Vh_.find(_vertex)->second;
        else return iGameVertexHandle(-1);
    }
    const iGameEdgeHandle Mesh::edgehandle(iGameEdge& _edge) const {
        if (Edge2Eh_.find(_edge) != Edge2Eh_.end()) return Edge2Eh_.find(_edge)->second;
        else return iGameEdgeHandle(-1);
    }
    const iGameFaceHandle Mesh::facehandle(iGameFace& _face) const {
        if (Face2Fh_.find(_face) != Face2Fh_.end()) return Face2Fh_.find(_face)->second;
        else return iGameFaceHandle(-1);
    }

    /*======================得到邻接关系============================*/
    // 顶点的邻接点
    // 先找邻接边，再找邻接点
    std::unordered_set<iGameVertexHandle> Mesh::NeighborVh(iGameVertexHandle _vh) {
        std::unordered_set<iGameVertexHandle> neighborvh;
        auto neighboreh = NeighborEh(_vh);
        // 存在邻接边是前提
        if (neighboreh.size()) {
            for (iGameEdgeHandle eh : neighboreh) {
                if (edges_[eh].vh1() != _vh) neighborvh.insert(edges_[eh].vh1());
                if (edges_[eh].vh2() != _vh) neighborvh.insert(edges_[eh].vh2());
            }
        }
        return neighborvh;
    }
    // 顶点的邻接边
    std::unordered_set<iGameEdgeHandle>& Mesh::NeighborEh(iGameVertexHandle _vh) {
        if (NeighborEhOfVertex_.count(_vh)) return NeighborEhOfVertex_[_vh];
        else return empty_ehs;               // 返回一个空的集合
    }
    // 顶点的邻接面
    std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameVertexHandle _vh) {
        if (NeighborFhOfVertex_.count(_vh)) return NeighborFhOfVertex_[_vh];
        else return empty_fhs;               // 返回一个空的集合
    }
    // 根据顶点的一条边得到对边的点： ·—·
    MeshKernel::iGameVertexHandle Mesh::NeighborVhFromEdge(iGameVertexHandle _vh, iGameEdgeHandle _eh) {
        assert(NeighborEh(_vh).count(_eh));     // 确保是邻接边
        iGameVertexHandle vh = (edges_[_eh].vh(0) == _vh) ? edges_[_eh].vh(1) : edges_[_eh].vh(0);
        return vh;
    }
    // 根据顶点的一条边得到对点的边： —·— ，四边形专用
    MeshKernel::iGameEdgeHandle Mesh::NeighborEhFromVertex(iGameEdgeHandle _eh, iGameVertexHandle _vh){
        assert(NeighborEh(_vh).count(_eh));     // 确保是邻接边
        if(!NeighborEh(_vh).count(_eh)) return EH{-1};
        if(NeighborEh(_vh).size() != 4) return EH{-1};// 奇异点直接返回-1
        std::vector<iGameEdgeHandle> n_ehs;
        auto n_fhs = NeighborFh(_eh);
        for(auto n_eh : NeighborEh(_vh)){
            if(n_eh != _eh) n_ehs.emplace_back(n_eh);
        }
        for(auto n_eh : n_ehs){
            bool not_neighbor = true;
            for(auto n_fh : NeighborFh(n_eh)){
                if(n_fhs.count(n_fh)){
                    not_neighbor = false;
                    break;
                }
            }
            if(not_neighbor) return n_eh;
        }
        // 如果找不到就返回-1
        return EH{-1};
    }

    // 边的邻接边
    // 两个顶点的所有邻接边去除当前边
    std::unordered_set<iGameEdgeHandle> Mesh::NeighborEh(iGameEdgeHandle _eh) {
        assert(edges_.count(_eh));                     // 保证该边handle存在
        std::unordered_set<iGameEdgeHandle> neighboreh;     // 保存输出的结果
        int k = 0;                                     // 遍历两个顶点
        while (k < 2) {
            iGameVertexHandle vh = edges_[_eh].vh(k);
            auto vhneighboreh = NeighborEh(vh);        // 得到点的邻接边
            for (iGameEdgeHandle eh : vhneighboreh) {
                if (eh != _eh) neighboreh.insert(eh);
            }
            ++k;
        }

        return neighboreh;
    }
    // 边的邻接面
    std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameEdgeHandle _eh) {
        if (NeighborFhOfEdge_.count(_eh)) return NeighborFhOfEdge_[_eh];
        else return empty_fhs;               // 返回一个空的集合
    }

    // 面的邻接面
    // 邻接面：有一条相同边
    std::unordered_set<iGameFaceHandle> Mesh::NeighborFh(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));                     // 保证该边handle存在
        std::unordered_set<iGameFaceHandle> neigborface;
        int k = 0;                                     // 遍历两个顶点
        size_t facesize = faces_[_fh].size();
        while (k < facesize) {
            iGameEdgeHandle eh = faces_[_fh].eh(k);
            auto ehneighborfh = NeighborFh(eh);        // 得到边的邻接面
            for (iGameFaceHandle fh : ehneighborfh) {
                if (fh != _fh) neigborface.insert(fh);
            }
            ++k;
        }
        return neigborface;
    }
    // 邻接面2：有一个公共顶点
    std::unordered_set<iGameFaceHandle> Mesh::Neighbor2Fh(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));                     // 保证该边handle存在
        std::unordered_set<iGameFaceHandle> neigborface;
        auto v_indices = faces_[_fh].getVertexHandle();
        for (auto& v_idx : v_indices) {
            auto adjF = NeighborFh(v_idx);
            for (iGameFaceHandle fh : adjF) {
                if (fh != _fh) neigborface.insert(fh);
            }
        }
        return neigborface;
    }

    // 根据面的一条边得到对边的面： 口|口 s
    iGameFaceHandle Mesh::NeighborFhFromEdge(iGameFaceHandle _fh, iGameEdgeHandle _eh){
        assert(NeighborFh(_eh).count(_fh));     // 保证是邻接边
        if (!NeighborFh(_eh).count(_fh)) return iGameFaceHandle{-1};
        for (iGameFaceHandle fh : NeighborFh(_eh)) {
            if (fh != _fh) return fh;
        }
    }

    // 根据面的一条边得到对边： |口| ，四边形专用，三角形或多边形会返回随机的另一条边
    iGameEdgeHandle Mesh::OppositeEhFromEdge(iGameFaceHandle _fh, iGameEdgeHandle _eh){
        assert(NeighborFh(_eh).count(_fh));     // 保证是邻接边
        auto& f = faces(_fh);
        if(f.size() == 4){
            for(auto n_eh : f.getEdgeHandle()){
                if(n_eh == _eh) continue;
                bool notNeighbor = true;
                for(auto n_n_eh : NeighborEh(_eh)){
                    if(n_n_eh == n_eh){
                        notNeighbor = false;
                        break;
                    }
                }
                if(notNeighbor) return n_eh;
            }
        }
        else{
            for(auto n_eh : f.getEdgeHandle()){
                if(n_eh == _eh) continue;
                return n_eh;
            }
        }
        return iGameEdgeHandle{-1};
    }

    // 根据两个顶点找到边，没找到就返回-1
    iGameEdgeHandle Mesh::getEhFromTwoVh(iGameVertexHandle vh1, iGameVertexHandle vh2) {
        iGameEdgeHandle tempEh{-1};
        for (auto eh : NeighborEh(vh1)) {
            auto& e = edges(eh);
            if (vh2 == (e.vh1() == vh1 ? e.vh2() : e.vh1())) {
                return eh;
            }
        }
        return tempEh;
    }

    iGameVertex Mesh::getEdgeVector(iGameEdge _e) {
        return vertices(_e.vh2()) - vertices(_e.vh1());
    }
}



// SurfaceMesh 定义
namespace MeshKernel {
    void SurfaceMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
                               const std::vector<std::vector<iGameVertexHandle>>& _elements) {
        std::vector<iGameVertexHandle>new_handle(_vertices.size());
        for (int i=0;i<_vertices.size();i++) {
            auto v =  _vertices[i];
            auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
            new_handle[i]=vh;
        }
        for (auto f : _elements) {
            //auto vh0 = AddVertex(iGameVertex(_vertices[f[0]].x(), v.y(), v.z()));
            if(f.size()==3)
            AddFace({new_handle[f[0]],new_handle[f[1]],new_handle[f[2]]});
            else if(f.size()==4)
            AddFace({new_handle[f[0]],new_handle[f[1]],new_handle[f[2]],new_handle[f[3]]});
        }
//        for (auto v : _vertices) {
//            auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
//
//        }
//        for (auto f : _elements) {
//            AddFace(f);
//        }
    }
    SurfaceMesh& SurfaceMesh::operator=(const SurfaceMesh& _surfacemesh) {
        if (this != &_surfacemesh) {
            Mesh::operator=(_surfacemesh);
        }
        return *this;
    }

    bool SurfaceMesh::isOnBoundary(iGameEdgeHandle eh) {
        auto fcnt = NeighborFh(eh).size();
        return fcnt == 1;
    }

    bool SurfaceMesh::isOnBoundary(iGameVertexHandle vh) {
        for (auto eh : NeighborEh(vh)) {
            if (isOnBoundary(eh)) {
                return true;
            }
        }
        return false;
    }

    bool SurfaceMesh::isOnBoundary(iGameFaceHandle fh) {
        auto face = faces(fh);
        for (auto eh : face.getEdgeHandle()) {
            if (isOnBoundary(eh)) {
                return true;
            }
        }
        return false;
    }

    bool SurfaceMesh::isTriangleMesh() {
        for (auto& fp : faces_) {
            if (fp.second.getVertexHandle().size() != 3) return false;
        }
        return true;
    }

    std::unordered_map<int, int> SurfaceMesh::updateAllHandles() {
        int vcnt = VertexSize(), fcnt = FaceSize();
        std::vector<MeshKernel::iGameVertex> newVertices;
        std::vector<std::vector<MeshKernel::iGameVertexHandle>> newFaces;
        std::unordered_map<int, int> mp;// old id to new id
        int idx = 0;
        for (auto& fp : allfaces()) {
            auto vhs = fp.second.getVertexHandle();
            for (auto& vh : vhs) {
                if (!mp.count(vh)) {
                    mp[vh] = idx++;
                    newVertices.push_back(vertices_[vh]);
                }
                vh = iGameVertexHandle(mp[vh]);
            }
            newFaces.push_back(vhs);
        }
        *this = MeshKernel::SurfaceMesh(newVertices, newFaces);
        return mp;
    }


    bool SurfaceMesh::isConnected(iGameVertexHandle vh1, iGameVertexHandle vh2) {
        for (auto vh : NeighborVh(vh1)) {
            if (vh == vh2) return true;
        }
        return false;
    }

    bool SurfaceMesh::isConnected(iGameEdgeHandle eh1, iGameEdgeHandle eh2) {
        if (!isValid(eh1) || !isValid(eh2)) return false;
        auto e1 = edges_[eh1];
        auto e2 = edges_[eh2];
        auto vh1 = e1.vh1(), vh2 = e1.vh2();
        auto vh3 = e2.vh1(), vh4 = e2.vh2();
        return (vh1 == vh3 || vh1 == vh4 || vh2 == vh3 || vh2 == vh4);
    }

    void SurfaceMesh::genNormal(iGameFaceHandle fh) {
        auto& face = faces(fh);
        auto vex = face.getVertexHandle();
        int n = vex.size();
        iGameVertex N;
        for (int i = 2; i < n; ++i) {
            auto& v0 = vertices(vex[0]);
            auto& v1 = vertices(vex[1]);
            auto& v2 = vertices(vex[2]);
            N += (v1 - v0).normalize() % (v2 - v0).normalize();
        }
        N.normalize();
        face.setNormal(N.x(), N.y(), N.z());
    }

    void SurfaceMesh::genNormal(iGameVertexHandle vh) {
        auto& v = vertices(vh);
        Eigen::Vector3d N = Eigen::Vector3d::Zero();
        auto adjFH = NeighborFh(vh);
        for (auto& fh : adjFH) {
            auto& face = faces(fh);
            Eigen::Vector3d faceN = Eigen::Vector3d(face.getNormalX(), face.getNormalY(), face.getNormalZ());
            N += faceN;
        }
        N /= adjFH.size();
        N.normalize();
        v.setNormal(N[0], N[1], N[2]);
    }

    void SurfaceMesh::genAllFacesNormal() {
        for (auto& fp : this->faces_) {
            genNormal(fp.first);
        }
    }

    void SurfaceMesh::genAllVerticesNormal() {
        this->genAllFacesNormal();
        for (auto& vp : this->vertices_) {
            genNormal(vp.first);
        }
    }

    void SurfaceMesh::genAllEdgesLength() {
        for (auto& ep : this->edges_) {
            genLength(ep.first);
        }
    }

    void SurfaceMesh::genLength(iGameEdgeHandle eh) {
        auto& e = edges(eh);
        auto vh1 = e.vh1();
        auto vh2 = e.vh2();
        auto& v1 = vertices(vh1);
        auto& v2 = vertices(vh2);
        Eigen::Vector3d vec(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
        e.setLength(vec.norm());
    }

    iGameEdgeHandle SurfaceMesh::getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2) {
        iGameEdgeHandle ret(-1);
        int vh_sum = vh1 + vh2;
        if (!isValid(vh1) || !isValid(vh2)) return ret;
        for (auto& eh : NeighborEh(vh1)) {
            auto& e = edges_[eh];
            if (vh_sum == e.vh1() + e.vh2()) {
                ret = eh;
                break;
            }
        }
        return ret;
    }

    size_t SurfaceMesh::getBoundaryVerticesCount() {
        size_t cnt = 0;
        for (auto& vp : vertices_) {
            if (isOnBoundary(vp.first))
                cnt++;
        }
        return cnt;
    }
    bool SurfaceMesh::isClosure() {

        for(auto i : this->alledges()){
            if(NeighborFh(i.first).size() == 1)
                return false;
        }
        return true;
    }

    double SurfaceMesh::getLength(iGameEdgeHandle eh) {
        assert(edges_.count(eh));
        auto& edge = edges_[eh];
        auto& v1 = vertices_[edge.vh1()];
        auto& v2 = vertices_[edge.vh2()];
        return (v1 - v2).norm();
    }

    double SurfaceMesh::getFaceArea(iGameFaceHandle fh){
        const auto& face = faces_[fh];
        auto vhs = face.getVertexHandle();
        if(vhs.size() == 3){
            double tri_area = 0;
            const auto& v1 = vertices_[vhs[0]];
            const auto& v2 = vertices_[vhs[1]];
            const auto& v3 = vertices_[vhs[2]];
            auto vec13 = v3 - v1;
            auto vec12 = v2 - v1;
            tri_area += (vec12 % vec13).norm();
            return tri_area/2;
        }
        else if(vhs.size() == 4){
            double quad_area = 0;
            const auto& v1 = vertices_[vhs[0]];
            const auto& v2 = vertices_[vhs[1]];
            const auto& v3 = vertices_[vhs[2]];
            const auto& v4 = vertices_[vhs[3]];
            auto vec13 = v3 - v1;
            auto vec12 = v2 - v1;
            auto vec14 = v4 - v1;
            quad_area += (vec12 % vec13).norm();
            quad_area += (vec14 % vec13).norm();
            return quad_area / 2;
        }
        else return 0;
    }

    void SurfaceMesh::addNoise(double level = 0.1) {
        const int I_MOD = 1E+6F;
        const double D_MOD = I_MOD;
        this->genAllVerticesNormal();
        this->genAllEdgesLength();
        int vcnt = this->vertices_.size();
        srand((unsigned)time(NULL));
        for (auto vp : vertices_) {
            auto vh = vp.first;
            auto& v = this->vertices(vh);
            double len = 0.f;
            size_t ecnt = 0;
            for (auto eh : NeighborEhOfVertex_[vh]) {
                len += (edges_[eh]).getLength();
                ecnt++;
            }
            if (ecnt == 0) continue;
            len /= ecnt;
            len *= level;
            double randX = (std::rand() % I_MOD) / D_MOD - 0.5f;
            double randY = (std::rand() % I_MOD) / D_MOD - 0.5f;
            double randZ = (std::rand() % I_MOD) / D_MOD - 0.5f;
            Eigen::Vector3d dir(v.getNormalX() + randX, v.getNormalY() + randY, v.getNormalZ() + randZ);
            dir.normalize();
            int flag = 2 * (std::rand() % 2) - 1;
            dir *= flag * len;
            v.setPosition(v.x() + dir.x(), v.y() + dir.y(), v.z() + dir.z());
        }
    }
}


// VolumeMesh 成员函数定义
namespace MeshKernel {
    void VolumeMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
                              const std::vector<std::vector<iGameVertexHandle>>& _elements) {
        for (auto v : _vertices) {
            auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
        }
        for (auto c : _elements) {
            AddCell(c);
        }
    }
    VolumeMesh& VolumeMesh::operator=(const VolumeMesh& _volumemesh) {
        if (this != &_volumemesh) {
            Mesh::operator=(_volumemesh);
        }
        this->empty_chs = _volumemesh.empty_chs;
        this->surface_faces_= _volumemesh.surface_faces_;

        this->cells_= _volumemesh.cells_;
        this->Cell2Ch_= _volumemesh.Cell2Ch_;

        this->NeighborChOfVertex_= _volumemesh.NeighborChOfVertex_;          //点的邻接体
        this-> NeighborChOfEdge_= _volumemesh.NeighborChOfEdge_;              //边的邻接体
        this-> NeighborChOfFace_= _volumemesh.NeighborChOfFace_;              //面的邻接体

        return *this;
    }
    iGameCell& VolumeMesh::cells(iGameCellHandle _ch) {
        assert(cells_.count(_ch));
        return cells_[_ch];
    }
    const iGameCell VolumeMesh::cells(iGameCellHandle _ch) const {
        assert(cells_.count(_ch));
        return cells_.find(_ch)->second;
    }
    const iGameCellHandle VolumeMesh::cellhandle(iGameCell& _cell) const {
        if (Cell2Ch_.find(_cell) != Cell2Ch_.end()) return Cell2Ch_.find(_cell)->second;
        else return iGameCellHandle(-1);
    }
    std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameVertexHandle _vh) {
        if (NeighborChOfVertex_.count(_vh)) return NeighborChOfVertex_[_vh];
        else return std::unordered_set<iGameCellHandle>();
    }
    std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameEdgeHandle _eh) {
        if (NeighborChOfEdge_.count(_eh)) return NeighborChOfEdge_[_eh];
        else return std::unordered_set<iGameCellHandle>();
    }
    std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameFaceHandle _fh) {
        if (NeighborChOfFace_.count(_fh)) return NeighborChOfFace_[_fh];
        else return std::unordered_set<iGameCellHandle>();
    }
    std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameCellHandle _ch) {
        assert(cells_.count(_ch));                     // 保证该边handle存在
        std::unordered_set<iGameCellHandle> neigborcell;
        int k = 0;                                     // 遍历两个顶点
        size_t facesize = cells_[_ch].faces_size();
        while (k < facesize) {
            iGameFaceHandle fh = cells_[_ch].fh(k);
            auto fhneighborch = NeighborCh(fh);        // 得到点的邻接边
            for (iGameCellHandle ch : fhneighborch) {
                if (ch != _ch) neigborcell.insert(ch);
            }
            ++k;
        }
        return neigborcell;
    }

    iGameCellHandle VolumeMesh::AddCell(const std::vector<iGameVertexHandle>& _vhs) {
        //std::cout << "加入的体中的各个点为 : " << std::endl;
        //for (int i = 0; i < _vhs.size(); i++) std::cout << _vhs[i] << " ";
        //std::cout << std::endl;
        int facecnt = _vhs.size() == 8 ? 6 : 4;
        int edgecnt = _vhs.size() == 8 ? 12 : 6;
        std::vector<iGameFaceHandle> fhs(facecnt, (iGameFaceHandle)0);
        std::vector<std::vector<int>> faceform(facecnt);
        if (facecnt == 6) {
            faceform = { {0,3,2,1},{0,4,7,3},{0,1,5,4},{4,5,6,7},{1,2,6,5},{2,3,7,6} };
        }
        else {
            faceform = { {0,1,3},{1,2,3},{2,0,3},{0,2,1} };
        }
        for (int i = 0; i < facecnt; ++i) {
            std::vector<iGameVertexHandle> facevertices(faceform[i].size());
            for (int j = 0; j < faceform[i].size(); ++j) {
                facevertices[j] = _vhs[faceform[i][j]];
            }
            fhs[i] = AddFace(facevertices);
        }
        std::vector<iGameEdgeHandle> ehs(edgecnt);
        ////////////////////// Test Begin
        //std::cout << "该体的每一个面有的边的数量 以及 fhs[] 中的数  : " << std::endl;
        //for (int i = 0; i < facecnt; i++) {
        //	std::cout << std::to_string(i) << " : " << fhs[i] << " ";
        //	std::cout<< faces_[fhs[i]].getEN() << std::endl;
        //}
        //std::cout << std::endl;
        ////////////////////// Test End
        if (edgecnt == 12) {
            ehs = { faces_[fhs[0]].eh(0),faces_[fhs[0]].eh(1),faces_[fhs[0]].eh(2),faces_[fhs[0]].eh(3),
                    faces_[fhs[3]].eh(0),faces_[fhs[3]].eh(1),faces_[fhs[3]].eh(2),faces_[fhs[3]].eh(3),
                    faces_[fhs[1]].eh(1), faces_[fhs[1]].eh(3), faces_[fhs[4]].eh(2), faces_[fhs[4]].eh(0) };
        }
        else {
            ehs = { faces_[fhs[3]].eh(0),faces_[fhs[3]].eh(1),faces_[fhs[3]].eh(2),
                    faces_[fhs[0]].eh(0), faces_[fhs[1]].eh(0), faces_[fhs[2]].eh(0) };
        }
        iGameCell c(_vhs, ehs, fhs);
        // 如果该体已经存在，则返回面的handle
        if (Cell2Ch_.count(c)) {
            //std::cout<<"该体已经存在 ." << std::endl;
            return Cell2Ch_[c];
        }
        // 否则
        else {
            iGameCellHandle ch = GenCellHandle();						// 生成一个新的handle
            cells_[ch] = c;											// 建立handle与该面之间的映射
            Cell2Ch_[c] = ch;										// 建立该面与handle之间的映射
            AddCell2Neighbor(ch);									// 将该面添加至面所包含的点和边的相邻面中
            return Cell2Ch_[c];										// 返回该面的handle
        }
    }

//    iGameCellHandle VolumeMesh::AddCell(const std::vector< std::vector<iGameVertexHandle>>& _vhs) {
//        if (_vhs.size() < 4) return iGameCellHandle(-1);
//        std::vector<MeshKernel::iGameVertexHandle> vertices;
//        std::vector<MeshKernel::iGameEdgeHandle> edges;
//        std::vector<MeshKernel::iGameFaceHandle> faces;
//        for (int i = 0; i < _vhs.size(); i++) {
//            faces.push_back(AddFace(_vhs[i]));
//            for (int j = 0; j < _vhs[i].size(); j++) {
//                if (j == 0) {
//                    auto eh = AddEdge(_vhs[i][_vhs[i].size() - 1], _vhs[i][j]);
//                    auto it = std::find(edges.begin(), edges.end(), eh);
//                    // 防止handle重复
//                    if (it == edges.end())
//                        edges.push_back(eh);
//                }
//                else {
//                    auto eh = AddEdge(_vhs[i][j], _vhs[i][j - 1]);
//                    auto it = std::find(edges.begin(), edges.end(), eh);
//                    // 防止handle重复
//                    if (it == edges.end())
//                        edges.push_back(eh);
//                }
//                auto it = std::find(vertices.begin(), vertices.end(), _vhs[i][j]);
//                // 防止handle重复
//                if (it == vertices.end())
//                    vertices.push_back(_vhs[i][j]);
//            }
//        }
//        iGameCell c(vertices, edges, faces);
//        iGameCellHandle ch = GenCellHandle();
//        cells_[ch] = c;
//        return ch;
//    }

    iGameVertexHandle VolumeMesh::DeleteVertex(const iGameVertexHandle& _vh) {
        if (!vertices_.count(_vh)) return iGameVertexHandle(-1);        // 如果该点不存在，则返回-1
        else {
            // 删除相邻边元素（删除边时自然删除了相邻面）
            auto ve = NeighborEhOfVertex_[_vh];
            for (iGameEdgeHandle eh : ve) {
                DeleteEdge(eh);
            }
            // 删除顶点元素和以其为根据的所有邻接关系
            Vertex2Vh_.erase(vertices_[_vh]);
            vertices_.erase(_vh);
            NeighborEhOfVertex_.erase(_vh);
            NeighborFhOfVertex_.erase(_vh);
            NeighborChOfVertex_.erase(_vh);
            return _vh;
        }
    }
    iGameEdgeHandle VolumeMesh::DeleteEdge(const iGameEdgeHandle& _eh) {
        if (!edges_.count(_eh)) return iGameEdgeHandle(-1);             // 如果该边不存在，则返回-1
        else {
            // 删除邻接关系
            iGameEdge e(edges_[_eh]);
            for (int i = 0; i < 2; ++i) {
                iGameVertexHandle ev = e.vh(i);
                NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
            }
            // 删除相邻面元素
            auto ef = NeighborFhOfEdge_[_eh];
            for (iGameFaceHandle fh : ef) {
                DeleteFace(fh);
            }
            // 删除边元素和以其为根据的所有邻接关系
            Edge2Eh_.erase(edges_[_eh]);
            edges_.erase(_eh);
            NeighborFhOfEdge_.erase(_eh);
            NeighborChOfEdge_.erase(_eh);
            return _eh;
        }
    }
    iGameFaceHandle VolumeMesh::DeleteFace(const iGameFaceHandle& _fh) {
        if (!faces_.count(_fh)) return iGameFaceHandle(-1);             // 如果该面不存在，则返回-1
        else {                                                     // 如果该面存在，则返回删除的这个面的handle
            // 删除邻接关系
            iGameFace f(faces_[_fh]);
            for (int i = 0; i < f.size(); ++i) {
                iGameVertexHandle fv = f.vh(i);
                iGameEdgeHandle fe = f.eh(i);
                NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
                NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
            }
            auto fc = NeighborChOfFace_[_fh];
            for (iGameCellHandle ch : fc) {
                DeleteCell(ch);
            }
            // 删除面元素
            Face2Fh_.erase(faces_[_fh]);
            faces_.erase(_fh);
            NeighborChOfFace_.erase(_fh);
            return _fh;
        }
    }
    iGameCellHandle VolumeMesh::DeleteCell(const iGameCellHandle& _ch) {
        if (!cells_.count(_ch)) return iGameCellHandle(-1);
        else {
            // 删除邻接关系
            iGameCell c(cells_[_ch]);
            std::vector<int> vsize = { 4,8 };
            std::vector<int> esize = { 6,12 };
            std::vector<int> fsize = { 4,6 };
            int volumeType = c.vertices_size() == 4 ? 0 : 1;
            for (int i = 0; i < vsize[volumeType]; ++i) {
                iGameVertexHandle cv = c.vh(i);
                NeighborChOfVertex_[cv].erase(_ch);
            }
            for (int i = 0; i < esize[volumeType]; ++i) {
                iGameEdgeHandle ce = c.eh(i);
                NeighborChOfEdge_[ce].erase(_ch);
            }
            for (int i = 0; i < fsize[volumeType]; ++i) {
                iGameFaceHandle cf = c.fh(i);
                NeighborChOfFace_[cf].erase(_ch);
            }
            // 删除面元素
            Cell2Ch_.erase(cells_[_ch]);
            cells_.erase(_ch);
            return _ch;
        }
    }

    void VolumeMesh::AddCell2Neighbor(const iGameCellHandle& _ch)
    {
        iGameCell c(cells_[_ch]);
        std::vector<int> vsize = { 4,8 };
        std::vector<int> esize = { 6,12 };
        std::vector<int> fsize = { 4,6 };
        int volumeType = c.vertices_size() == 4 ? 0 : 1;
        for (int i = 0; i < vsize[volumeType]; ++i) {
            NeighborChOfVertex_[c.vh(i)].insert(_ch);
        }
        for (int i = 0; i < esize[volumeType]; ++i) {
            NeighborChOfEdge_[c.eh(i)].insert(_ch);
        }
        for (int i = 0; i < fsize[volumeType]; ++i) {
            NeighborChOfFace_[c.fh(i)].insert(_ch);
        }
    }

    bool VolumeMesh::isConnected(iGameFaceHandle& fh1, iGameFaceHandle& fh2) {
        auto& face1 = faces_[fh1];
        auto& face2 = faces_[fh2];
        const auto& ehs1 = face1.getEdgeHandle();
        for (const auto& eh : face2.getEdgeHandle()) {
            if (std::find(ehs1.begin(), ehs1.end(), eh) != ehs1.end()) {
                return true;
            }
        }
        return false;
    }

    bool VolumeMesh::isConnected(iGameEdgeHandle& eh1, iGameEdgeHandle& eh2) {// eh1 == eh2 is not connected
        if (eh1 == eh2) return false;
        return edges_[eh1].vh1() == edges_[eh2].vh1() || edges_[eh1].vh1() == edges_[eh2].vh2() ||
               edges_[eh1].vh2() == edges_[eh2].vh1() || edges_[eh1].vh2() == edges_[eh2].vh2();
    }

    bool VolumeMesh::isConnected(iGameVertexHandle& vh1, iGameVertexHandle& vh2) {
        const auto& adjvhs = NeighborVh(vh1);
        for (const auto& adjvh : adjvhs) {
            if (adjvh == vh2) return true;
        }
        return false;
    }
    bool VolumeMesh::isConnected(const iGameVertexHandle& vh1, const iGameVertexHandle& vh2) {
        const auto& adjvhs = NeighborVh(vh1);
        for (const auto& adjvh : adjvhs) {
            if (adjvh == vh2) return true;
        }
        return false;
    }

    bool VolumeMesh::isOnBoundary(iGameCellHandle ch) {
        return NeighborCh(ch).size() < 6;
    }

    bool VolumeMesh::isOnBoundary(iGameFaceHandle fh) {
        return NeighborCh(fh).size() == 1;// 0 isnot on the boundary
    }

    bool VolumeMesh::isOnBoundary(iGameEdgeHandle eh) {
        return NeighborCh(eh).size() < 4;
    }

    bool VolumeMesh::isOnBoundary(iGameVertexHandle vh) {
        return NeighborCh(vh).size() < 8;
    }

    iGameVertex VolumeMesh::getCenter(iGameCellHandle ch) {
        iGameVertex center(0, 0, 0);
        const auto& cell = cells_[ch];
        auto vhs = cell.getVertexHandle();
        for (auto& vh : vhs) {
            center = center + vertices_[vh];
        }
        return center / vhs.size();
    }

    double VolumeMesh::getQuadArea(iGameFaceHandle fh) {
        const auto& face = faces_[fh];
        auto vhs = face.getVertexHandle();
        assert(vhs.size() == 4);
        double quad_area = 0;
        const auto& v1 = vertices_[vhs[0]];
        const auto& v2 = vertices_[vhs[1]];
        const auto& v3 = vertices_[vhs[2]];
        const auto& v4 = vertices_[vhs[3]];
        auto vec13 = v3 - v1;
        auto vec12 = v2 - v1;
        auto vec14 = v4 - v1;
        quad_area += (vec12 % vec13).norm() * 0.5;
        quad_area += (vec14 % vec13).norm() * 0.5;
        quad_area += (vec12 % vec14).norm() * 0.5;
        auto vec23 = v3 - v2;
        auto vec24 = v4 - v2;
        quad_area += (vec23 % vec24).norm() * 0.5;
        return quad_area / 4;
    }

    iGameVertex VolumeMesh::getQuadNormal(iGameFaceHandle fh) {// 返回的是经过面积加权的
        iGameVertex N(0, 0, 0);
        const auto& face = faces_[fh];
        auto vhs = face.getVertexHandle();
        auto chs = NeighborCh(fh);
        if (chs.size() != 1) return N;
        assert(vhs.size() == 4);
        auto cell_cenetr = getCenter(*chs.begin());
        //double quad_area = 0;
        const auto& v1 = vertices_[vhs[0]];
        const auto& v2 = vertices_[vhs[1]];
        const auto& v3 = vertices_[vhs[2]];
        const auto& v4 = vertices_[vhs[3]];
        auto face_center = (v1 + v2 + v3 + v4) / 4;
        auto out_dir = face_center - cell_cenetr;
        out_dir = out_dir.normalize();
        auto vec13 = v3 - v1;
        auto vec12 = v2 - v1;
        auto vec14 = v4 - v1;
        auto vec23 = v3 - v2;
        auto vec24 = v4 - v2;

        auto N123 = vec12 % vec13;
        double area1 = N123.norm() * 0.5;
        N123 = N123.normalize();
        if (N123 * out_dir < 0) N123 = N123 * -1;

        auto N134 = vec14 % vec13;
        double area2 = N134.norm() * 0.5;
        N134 = N134.normalize();
        if (N134 * out_dir < 0) N134 = N134 * -1;

        auto N124 = vec12 % vec14;
        double area3 = N124.norm() * 0.5;
        N124 = N124.normalize();
        if (N124 * out_dir < 0) N124 = N124 * -1;

        auto N234 = vec23 % vec24;
        double area4 = N234.norm() * 0.5;
        N234 = N234.normalize();
        if (N234 * out_dir < 0) N234 = N234 * -1;
        N = (N123 * area1 + N134 * area2 + N124 * area3 + N234 * area4);
        //N = N.normalize();
        return N;
    }

    double VolumeMesh::getLength(iGameEdgeHandle eh) {
        auto& v1 = vertices_[edges_[eh].vh1()];
        auto& v2 = vertices_[edges_[eh].vh2()];
        return (v1 - v2).norm();
    }

    void VolumeMesh::updateAllHandles() {
        std::vector<MeshKernel::iGameVertex> newVertices;
        std::vector<std::vector<MeshKernel::iGameVertexHandle>> newCells;
        std::unordered_map<int, int> mp;// old id to new id
        int idx = 0;
        for (auto& cp : allcells()) {
            auto vhs = cp.second.getVertexHandle();
            for (auto& vh : vhs) {
                if (!mp.count(vh)) {
                    mp[vh] = idx++;
                    newVertices.push_back(vertices_[vh]);
                }
                vh = iGameVertexHandle(mp[vh]);
            }
            newCells.push_back(vhs);
        }

        *this = MeshKernel::VolumeMesh(newVertices, newCells);

    }

    void VolumeMesh::genNormal(iGameFaceHandle fh) {
        auto& face = faces(fh);
        auto vex = face.getVertexHandle();
        std::vector<Eigen::Vector3d> pos(3, Eigen::Vector3d::Zero());
        for (int i = 0; i < 3; ++i) {
            auto& v = vertices(vex[i]);
            pos[i][0] = v.x();
            pos[i][1] = v.y();
            pos[i][2] = v.z();
        }
        Eigen::Vector3d N = (pos[1] - pos[0]).cross(pos[2] - pos[0]);
        N.normalize();
        face.setNormal(N[0], N[1], N[2]);
    }

    void VolumeMesh::genNormal(iGameVertexHandle vh) {
        auto& v = vertices(vh);
        Eigen::Vector3d N = Eigen::Vector3d::Zero();
        auto adjFH = NeighborFh(vh);
        for (auto& fh : adjFH) {
            auto& face = faces(fh);
            Eigen::Vector3d faceN = Eigen::Vector3d(face.getNormalX(), face.getNormalY(), face.getNormalZ());
            N += faceN;
        }
        N /= adjFH.size();
        N.normalize();
        v.setNormal(N[0], N[1], N[2]);
    }

    void VolumeMesh::genAllFacesNormal() {
        for (auto& fp : this->faces_) {
            genNormal(fp.first);
        }
    }

    void VolumeMesh::genAllVerticesNormal() {
        this->genAllFacesNormal();
        for (auto& vp : this->vertices_) {
            genNormal(vp.first);
        }
    }

}

namespace MeshKernel {

}

namespace MeshKernel {

}


