#pragma once
#include "Mesh.h"
#include "stl_parser.h"
#include "stl_write.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <array>
#include <map>
#include <float.h>

namespace MeshKernel {
	class IO {
	public:
		SurfaceMesh ReadObjFile(const std::string& _InputFile);
		SurfaceMesh ReadOffFile(const std::string& _InputFile);
		SurfaceMesh ReadStlFile(const std::string& _InputFile);
		SurfaceMesh ReadVtkFile(const std::string& _InputFile, std::vector<array<int, 2>>& quadconstrains);
		VolumeMesh ReadMeshFile(const std::string& _InputFile);
        VolumeMesh ReadVtkFile_Volume(const std::string& _InputFile);
        VolumeMesh m_read_mesh_file(const std::string& _InputFile);
		bool WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		bool WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		bool WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		bool WriteMeshFile(const VolumeMesh& _mesh, const std::string& _OutputFile);
		bool WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		std::string WriteOffString(const SurfaceMesh& _mesh);
		VolumeMesh ReadMeshFileFromStr(const std::string& data);
	private:
		void ReOrderVertexHandle(const SurfaceMesh& _mesh);
		std::vector<iGameVertexHandle> reorderedvh_;                        // 重排顶点
		std::unordered_map<iGameVertexHandle, std::size_t> newvh_;          // 新的顶点handle
	};
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadObjFile(const std::string& _InputFile) {
	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::vector<std::vector<double>> normals;
	std::vector<std::vector<double>> uvs;
	std::unordered_map<int, int> V2N;// vertex to normal
	std::unordered_map<int, int> V2T;// vertex to uv
	std::string line;

	std::cout << "Reading " << _InputFile << " File" << std::endl;
	// std::cout << inputfile.good() << std::endl;
	while (inputfile) {
		line.clear();
		getline(inputfile, line);
		if (line[0] == '#') {
			continue;// 注释
		}
		std::stringstream linestream;
		linestream.str(line);

		std::string flag;
		linestream >> flag;
		if (flag == "v") {
			double x, y, z;
			linestream >> x >> y >> z;
			vertices.push_back(iGameVertex(x, y, z));
		}
		else if (flag == "f") {
			// f 1575/1514/1569 1581/1520/1575 1576/1515/1570
			// f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
			std::vector<std::string> vex;
			std::string tmp;
			while (linestream >> tmp) vex.push_back(tmp);
			auto n = vex.size();
			//std::cout <<"n: "<<n << std::endl;
			// 创建面
			std::vector<iGameVertexHandle> face(n);
			for (size_t i = 0; i < n; i++) {
				size_t idx = 0;
				while (idx < vex[i].length() && isdigit(vex[i][idx])) idx++;
				int vh = std::stoi(vex[i].substr(0, idx)) - 1;// 注意 obj 文件 v 从1开始
				face[i] = (iGameVertexHandle)(vh);
				if (idx != vex[i].length()) {// 指定了纹理坐标
					// 注意：纹理坐标存在复用！法向量存在复用！顶点、纹理坐标与法向量个数均可互相不等
					size_t beg = idx;
                    int s = 0;// 斜杠个数
					while (beg < vex[i].length() && !isdigit(vex[i][beg])){
                        if(vex[i][beg] == '/') s++;
                        beg++;
                    }
                    size_t end = beg;
                    // 斜杠超过一个，说明纹理省略了，直接略过
                    if(s < 2){
                        while (end < vex[i].length() && isdigit(vex[i][end])) end++;
                    }
					if (!V2T.count(vh) && s < 2) {
						int uv_idx = std::stoi(vex[i].substr(beg, end - beg)) - 1;// 注意 obj 文件 vt 从1开始
						V2T[vh] = uv_idx;
                        //cout<<"stoi: "<<std::stoi(vex[i].substr(beg, end - beg))<<"  ";
                        //cout<<"uv_idx: "<<uv_idx<<endl;
					}
					if (end != vex[i].length() && !V2N.count(vh)) {// 指定了法向量
						beg = end;
						while (beg < vex[i].length() && !isdigit(vex[i][beg])) beg++;
						end = beg;
						while (end < vex[i].length() && isdigit(vex[i][end])) end++;
						int n_idx = std::stoi(vex[i].substr(beg, end - beg)) - 1;// 注意 obj 文件 vn 从1开始
						V2N[vh] = n_idx;
//                        cout<<"stoi: "<<std::stoi(vex[i].substr(beg, end - beg))<<endl;
//                        cout<<"n_idx: "<<n_idx<<endl;
					}
				}
			}
			faces.push_back(face);
		}
		else if (flag == "vt") {
			double u, v;
			linestream >> u >> v;
			uvs.push_back({ u, v });
		}
		else if (flag == "vn") {
			double x, y, z;
			linestream >> x >> y >> z;
			normals.push_back({ x, y, z });
		}
	}
	//printf("read file success, fcnt: %d, vcnt: %d, vtcnt: %d, vncnt: %d\n", faces.size(), vertices.size(), uvs.size(), normals.size());
	if (!normals.empty()) {
		int ncnt = normals.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int nidx = V2N[i];
            //cout<<"nidx: "<<nidx<<endl;
			//if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
			assert(nidx >= 0 && nidx < ncnt);
			vertices[i].setNormal(normals[nidx]);
		}
	}
	if (!uvs.empty()) {
		int uvcnt = uvs.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int uvidx = V2T[i];
			//if (uvidx < 0 || uvidx >= uvcnt) printf("error: uvidx = %d\n", uvidx);// debug 用
			assert(uvidx >= 0 && uvidx < uvcnt);
			//vertices[i].setUV(uvs[uvidx]);
		}
	}

	// 把模型位置移到[-1,1]范围内，中心点移到原点
	double minX = DBL_MAX; double maxX = -DBL_MAX;
	double minY = DBL_MAX; double maxY = -DBL_MAX;
	double minZ = DBL_MAX; double maxZ = -DBL_MAX;
	for (auto& v : vertices) {
		minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
		minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
		minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
	}
	double centerX = (minX + maxX) / 2;
	double centerY = (minY + maxY) / 2;
	double centerZ = (minZ + maxZ) / 2;
	double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
	for (auto& v : vertices) {
		v.x() -= centerX;
		v.y() -= centerY;
		v.z() -= centerZ;
		v = v / radius;
	}
	auto mesh = SurfaceMesh(vertices, faces);
	//if (!normals.empty()) mesh.setHasVN(true);
	//if (!uvs.empty()) mesh.setHasVT(true);
	inputfile.close();
	return mesh;
}


MeshKernel::SurfaceMesh MeshKernel::IO::ReadOffFile(const std::string& _InputFile) {
	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::vector<std::vector<double>> normals;
	std::vector<std::vector<double>> uvs;
	std::unordered_map<int, int> V2N;// vertex to normal
	std::unordered_map<int, int> V2T;// vertex to uv
	std::string line;
	int v_size, f_size, e_size;

	std::cout << "Reading " << _InputFile << " File" << std::endl;
	line.clear();
	getline(inputfile, line);
	std::stringstream linestream;

	if (line == "OFF") {
		line.clear();
		getline(inputfile, line);
		while(line.size() == 0 || line[0] == '#')
            getline(inputfile, line);
		linestream.str(line);
		linestream >> v_size >> f_size >> e_size;
	}
	for (int i = 0; i < v_size; i++) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		double x, y, z;
		linestream >> x >> y >> z;
		vertices.push_back(iGameVertex(x, y, z));
		//std::cout << x << " " << y << " " << z << std::endl;
	}
	for (int i = 0; i < f_size; i++) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		int v;
		linestream >> v;
		//std::cout << v << " " << x << " " << y << " " << z << std::endl;
		std::vector<iGameVertexHandle> face(v);
		for (int i = 0; i < v; i++) {
			int temp;
			linestream >> temp;
			face[i] = (iGameVertexHandle)(temp);
		}
		faces.push_back(face);
	}
	printf("read file success, fcnt: %ld, vcnt: %ld, vtcnt: %ld, vncnt: %ld\n", faces.size(), vertices.size(), uvs.size(), normals.size());
	if (!normals.empty()) {
		int ncnt = normals.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int nidx = V2N[i];
			//if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
			assert(nidx >= 0 && nidx < ncnt);
			vertices[i].setNormal(normals[nidx]);
		}
	}
	if (!uvs.empty()) {
		int uvcnt = uvs.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int uvidx = V2T[i];
			//if (uvidx < 0 || uvidx >= uvcnt) printf("error: uvidx = %d\n", uvidx);// debug 用
			assert(uvidx >= 0 && uvidx < uvcnt);
			//vertices[i].setUV(uvs[uvidx]);
		}
	}

	// 把模型位置移到[-1,1]范围内，中心点移到原点
    double minX = DBL_MAX; double maxX = -DBL_MAX;
    double minY = DBL_MAX; double maxY = -DBL_MAX;
    double minZ = DBL_MAX; double maxZ = -DBL_MAX;
	for (auto& v : vertices) {
		minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
		minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
		minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
	}
	double centerX = (minX + maxX) / 2;
	double centerY = (minY + maxY) / 2;
	double centerZ = (minZ + maxZ) / 2;
	double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
	for (auto& v : vertices) {
		v.x() -= centerX;
		v.y() -= centerY;
		v.z() -= centerZ;
		v = v / radius;
	}

	auto mesh = SurfaceMesh(vertices, faces);
	//if (!normals.empty()) mesh.setHasVN(true);
	//if (!uvs.empty()) mesh.setHasVT(true);
	inputfile.close();
	return mesh;
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadStlFile(const std::string& _InputFile) {

	char const* file_path = _InputFile.c_str();

	size_t file_size;

	FILE* file = fopen(file_path, "rb");
	if (!file) {
		fprintf(stderr, "file %s does not exist\n", file_path);
		exit(EXIT_FAILURE);
	}

	fseek(file, 0, SEEK_END);
	file_size = ftell(file);
	rewind(file);

	char* file_contents = (char*)malloc(file_size);
	if (fread(file_contents, 1, file_size, file) != file_size) {
		fprintf(stderr, "I/O error while reading %s\n", file_path);
		exit(EXIT_FAILURE);
	}
	fclose(file);

	size_t buf_size;
	if (one_stl_buf_size(&buf_size, file_contents, file_size, ONE_STL_NVVV)) {
		fprintf(stderr, "file \"%s\" is not valid stl data\n", file_path);
		exit(EXIT_FAILURE);
	}

	float* buf = (float*)malloc(buf_size);
	size_t trig_count = one_stl_parse(buf, file_contents, ONE_STL_NVVV);

	printf("triangle count: %ld\n", trig_count);

	typedef struct {
		float normal[3];
		float v0[3];
		float v1[3];
		float v2[3];
	} trig_nvvv_s;
	trig_nvvv_s* triangles = (trig_nvvv_s*)buf;
	float EPS = 0;
	//for (size_t i = 0; i < trig_count; ++i){
	//    EPS=std::max(EPS,(triangles + i)->v0[0]);
	//    EPS=std::max(EPS,(triangles + i)->v0[1]);
	//    EPS=std::max(EPS,(triangles + i)->v0[2]);
	//}
	//EPS/=100000;
	//EPS = 0;

	struct PointEps {
		PointEps() {}
		PointEps(float a, float b, float c, float EPS) {
			v0[0] = a; v0[1] = b;
			v0[2] = c;
			this->EPS = EPS;
		}
		float v0[3];
		float EPS;
		bool operator < (const PointEps& other) const {
			for (int i = 0; i < 3; i++) {
				if (abs(v0[i] - other.v0[i]) > EPS) {
					return  v0[i] + EPS < other.v0[i];
				}
			}
			return false;
		}
	};

	std::map<PointEps, int>mp;
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::function<int(float, float, float)> getHandleId = [&](float x, float y, float z) {
		if (mp.count(PointEps(x, y, z, EPS))) {
			return mp[PointEps(x, y, z, EPS)];
		}
		int cnt = mp.size();
		mp[PointEps(x, y, z, EPS)] = cnt;
		vertices.push_back(MeshKernel::iGameVertex{ x,y,z });
		return cnt;
	};
	std::vector<std::vector<float> >normal;
	for (size_t i = 0; i < trig_count; ++i) {
		faces.push_back({ (iGameVertexHandle)getHandleId((triangles + i)->v0[0],(triangles + i)->v0[1],(triangles + i)->v0[2]),
						(iGameVertexHandle)getHandleId((triangles + i)->v1[0],(triangles + i)->v1[1],(triangles + i)->v1[2]),
						(iGameVertexHandle)getHandleId((triangles + i)->v2[0],(triangles + i)->v2[1],(triangles + i)->v2[2])
			});
		normal.push_back(std::vector<float>{ (triangles + i)->normal[0], (triangles + i)->normal[1], (triangles + i)->normal[2] });
	}
	std::cout << vertices.size() << " " << faces.size() << std::endl;

	// 把模型位置移到[-1,1]范围内，中心点移到原点
    double minX = DBL_MAX; double maxX = -DBL_MAX;
    double minY = DBL_MAX; double maxY = -DBL_MAX;
    double minZ = DBL_MAX; double maxZ = -DBL_MAX;
	for (auto& v : vertices) {
		minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
		minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
		minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
	}
	double centerX = (minX + maxX) / 2;
	double centerY = (minY + maxY) / 2;
	double centerZ = (minZ + maxZ) / 2;
	double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
	for (auto& v : vertices) {
		v.x() -= centerX;
		v.y() -= centerY;
		v.z() -= centerZ;
		v = v / radius;
	}

	auto mesh = SurfaceMesh(vertices, faces);
	std::cout << trig_count << " " << faces.size() << std::endl;
	//for(int i=0;i<faces.size();i++){
	//    mesh.faces(FaceHandle(i)).setNormal(normal[i][0], normal[i][1], normal[i][2]);
	//}   
	free(file_contents);
	free(buf);
	return mesh;
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadVtkFile(const std::string& _InputFile, std::vector<array<int, 2>>& quadconstrains)
{
	std::vector<array<double, 3>> triverts;
	std::vector<array<uint32_t, 3>> triangles;
	FILE* fp = fopen(_InputFile.c_str(), "rb");

	char buffer[256], buffer2[256];
	std::map<int, std::map<int, std::string> > physicals[4];

	if (!fgets(buffer, sizeof(buffer), fp)) {
		fclose(fp);
	} // version line
	if (!fgets(buffer, sizeof(buffer), fp)) {
		fclose(fp);
	} // title

	if (fscanf(fp, "%s", buffer) != 1) // ASCII or BINARY
		fprintf(stdout, "Failed reading buffer\n");
	bool binary = false;
	if (!strcmp(buffer, "BINARY")) binary = true;

	if (fscanf(fp, "%s %s", buffer, buffer2) != 2) {
		fclose(fp);
	}

	bool unstructured = false;
	if (!strcmp(buffer, "DATASET") && !strcmp(buffer2, "UNSTRUCTURED_GRID"))
		unstructured = true;

	if ((strcmp(buffer, "DATASET") && strcmp(buffer2, "UNSTRUCTURED_GRID")) ||
		(strcmp(buffer, "DATASET") && strcmp(buffer2, "POLYDATA"))) {
		fprintf(stdout, "VTK reader can only read unstructured or polydata datasets\n");
		fclose(fp);
	}

	// read mesh vertices
	int numVertices;
	if (fscanf(fp, "%s %d %s\n", buffer, &numVertices, buffer2) != 3) ;
	if (strcmp(buffer, "POINTS") || !numVertices) {
		fprintf(stdout, "No points in dataset\n");
		fclose(fp);
	}
	int datasize;
	if (!strcmp(buffer2, "double"))
		datasize = sizeof(double);
	else if (!strcmp(buffer2, "float"))
		datasize = sizeof(float);
	else {
		fprintf(stdout, "VTK reader only accepts float or double datasets\n");
		fclose(fp);
	}
	fprintf(stdout, "Reading %d points\n", numVertices);
	std::vector<array<double, 3>> vertices(numVertices);
	for (int i = 0; i < numVertices; i++) {
		double xyz[3];
		if (binary) {
			if (datasize == sizeof(float)) {
				float f[3];
				if (fread(f, sizeof(float), 3, fp) != 3) {
					fclose(fp);
				}

				for (int j = 0; j < 3; j++) xyz[j] = f[j];
			}
			else {
				if (fread(xyz, sizeof(double), 3, fp) != 3) {
					fclose(fp);
				}

			}
		}
		else {
			if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
				fclose(fp);
			}
		}
		vertices[i] = array<double, 3>{xyz[0], xyz[1], xyz[2]};
	}

	// read mesh elements
	int numElements, totalNumInt;
	if (fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3) {
		fclose(fp);
	}

	bool haveCells = true;
	bool haveLines = false;
	if (!strcmp(buffer, "CELLS") && numElements > 0)
		fprintf(stdout, "Reading %d cells\n", numElements);
	else if (!strcmp(buffer, "POLYGONS") && numElements > 0)
		fprintf(stdout, "Reading %d polygons\n", numElements);
	else if (!strcmp(buffer, "LINES") && numElements > 0) {
		haveCells = false;
		haveLines = true;
		fprintf(stdout, "Reading %d lines\n", numElements);
	}
	else {
		fprintf(stdout, "No cells or polygons in dataset\n");
		fclose(fp);
	}


	if (haveCells) {
		std::vector<std::vector<int> > cells(numElements);
		for (std::size_t i = 0; i < cells.size(); i++) {
			int numVerts, n[100];
			if (binary) {
				if (fread(&numVerts, sizeof(int), 1, fp) != 1) {
					fclose(fp);
				}

				if ((int)fread(n, sizeof(int), numVerts, fp) != numVerts) {
					fclose(fp);
				}

			}
			else {
				if (fscanf(fp, "%d", &numVerts) != 1) {
					fclose(fp);
				}
				for (int j = 0; j < numVerts; j++) {
					if (fscanf(fp, "%d", &n[j]) != 1) {
						fclose(fp);
					}
				}
			}
			for (int j = 0; j < numVerts; j++) {
				if (n[j] >= 0 && n[j] < (int)vertices.size())
					cells[i].push_back(n[j]);
				else
					fprintf(stdout, "Wrong node index %d\n", n[j]);
			}
		}

		if (unstructured) {
			if (fscanf(fp, "%s %d\n", buffer, &numElements) != 2) {
				fclose(fp);

			}
			if (strcmp(buffer, "CELL_TYPES") || numElements != (int)cells.size()) {
				fprintf(stdout, "No or invalid number of cells types\n");
				fclose(fp);

			}
			for (std::size_t i = 0; i < cells.size(); i++) {
				int type;
				if (binary) {
					if (fread(&type, sizeof(int), 1, fp) != 1) {
						fclose(fp);

					}

				}
				else {
					if (fscanf(fp, "%d", &type) != 1) {
						fclose(fp);

					}
				}
				switch (type) {

					// first order elements
				case 3: {
					quadconstrains.push_back({ cells[i][0],cells[i][1] });
					break;
				}

				case 5: {
					triangles.push_back({ (uint32_t)cells[i][0],(uint32_t)cells[i][1],(uint32_t)cells[i][2] });
					break;
				}

				default: fprintf(stdout, "Unknown type of cell %d\n", type); break;
				}
			}
		}
		else {
			for (std::size_t i = 0; i < cells.size(); i++) {
				int nbNodes = (int)cells[i].size();
				switch (nbNodes) {

				case 2: {
					quadconstrains.push_back({ cells[i][0],cells[i][1] });
					break;
				}
				case 3: {
					triangles.push_back({ (uint32_t)cells[i][0],(uint32_t)cells[i][1],(uint32_t)cells[i][2] });
					break;
				}
				default:
					fprintf(stdout, "Unknown type of mesh element with %d nodes\n", nbNodes);
					break;
				}
			}
		}
	}
	else if (haveLines) {

	}

	triverts = vertices;

	fclose(fp);

	std::vector<iGameVertex> vertice;
	std::vector<std::vector<iGameVertexHandle>> faces;
	//std::vector<array<double, 3>> triverts;
	//std::vector<array<uint32_t, 3>> triangles;
	for (auto v : triverts) {
		vertice.emplace_back(iGameVertex{v[0],v[1],v[2]});
	}
	for (auto f : triangles) {
		faces.emplace_back(std::vector{iGameVertexHandle{(int)f[0]},iGameVertexHandle{(int)f[1]},iGameVertexHandle{(int)f[2]}});
	}

	auto mesh = SurfaceMesh(vertice, faces);
	return mesh;

}

std::string MeshKernel::IO::WriteOffString(const SurfaceMesh& _mesh) {
	std::stringstream ret;
	ReOrderVertexHandle(_mesh);
	ret << "OFF" << std::endl;
	ret << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << _mesh.alledges().size() << std::endl;
	for (iGameVertexHandle vh : reorderedvh_) {
		iGameVertex v(_mesh.vertices(vh));
		ret << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		ret << f.second.size();
		for (int i = 0; i < f.second.size(); ++i) {
			ret << " " << newvh_[f.second.vh(i)];
		}
		ret << std::endl;
	}
	return std::string(ret.str());
}

bool MeshKernel::IO::WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	std::ofstream outputfile(_OutputFile, std::ios::out);
//	ReOrderVertexHandle(_mesh);
//	for (iGameVertexHandle vh : reorderedvh_) {
//		iGameVertex v(_mesh.vertices(vh));
//		outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
//	}
//	auto allf = _mesh.allfaces();
//	for (auto f : allf) {
//		outputfile << "f";
//		for (int i = 0; i < f.second.size(); ++i) {
//			outputfile << " " << newvh_[f.second.vh(i)] + 1;
//		}
//		outputfile << std::endl;
//	}

    // TODO:只是为了布尔运算debug用，后面需要改回来！
    for(int i=0;i<_mesh.allvertices().size();i++){
        iGameVertex v(_mesh.vertices(iGameVertexHandle{i}));
        outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }
    for(int i=0;i<_mesh.allfaces().size();i++){
        auto& f = _mesh.faces(iGameFaceHandle{i});
        outputfile << "f";
        for (int j = 0; j < f.size(); ++j) {
            outputfile << " " << f.vh(j) + 1;
        }
        outputfile << std::endl;
    }

	outputfile.close();
	return true;
}

bool MeshKernel::IO::WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	std::ofstream outputfile(_OutputFile, std::ios::out);
	ReOrderVertexHandle(_mesh);
	outputfile << "OFF" << std::endl;
	outputfile << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << _mesh.alledges().size() << std::endl;
	for (iGameVertexHandle vh : reorderedvh_) {
		iGameVertex v(_mesh.vertices(vh));
		outputfile << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		outputfile << f.second.size();
		for (int i = 0; i < f.second.size(); ++i) {
			outputfile << " " << newvh_[f.second.vh(i)];
		}
		outputfile << std::endl;
	}
	outputfile.close();
	return true;
}

bool MeshKernel::IO::WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	STLWrite c;
	for (auto i : _mesh.allfaces()) {
		printf("size : %ld\n", i.second.getSortedVertexHandle().size());

		if (i.second.getSortedVertexHandle().size() == 3) {
			for (int j = 0; j < 3; j++) {
				auto point = _mesh.vertices(i.second.vh(j));
				c.add(point.x(), point.y(), point.z());
			}
		}
		else if (i.second.getSortedVertexHandle().size() == 4) {
			for (int j = 0; j < 3; j++) {
				auto point = _mesh.vertices(i.second.vh(j));
				c.add(point.x(), point.y(), point.z());
			}

			auto point = _mesh.vertices(i.second.vh(0));
			c.add(point.x(), point.y(), point.z());
			point = _mesh.vertices(i.second.vh(2));
			c.add(point.x(), point.y(), point.z());
			point = _mesh.vertices(i.second.vh(3));
			c.add(point.x(), point.y(), point.z());

		}
	}
	c.write(_OutputFile);
	return true;
}

std::vector<std::string> SplitFileName(const std::string& fileName)
{
	// JFR DO NOT CHANGE TO std::vector<std::string> s(3), it segfaults while
	// destructor si called
	std::vector<std::string> s;
	s.resize(3);
	if (fileName.size()) {
		// returns [path, baseName, extension]
		int idot = (int)fileName.find_last_of('.');
		int islash = (int)fileName.find_last_of("/\\");
		if (idot == (int)std::string::npos) idot = -1;
		if (islash == (int)std::string::npos) islash = -1;
		if (idot > 0) s[2] = fileName.substr(idot);
		if (islash > 0) s[0] = fileName.substr(0, islash + 1);
		s[1] =
			fileName.substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
	}
	return s;
}

bool MeshKernel::IO::WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	// 还不支持约束线的输出

	std::vector<std::array<double, 3>>verts;
	std::vector<std::vector<uint32_t>> faces;
	std::vector<std::array<int, 2>> constrains;

	for (int i = 0; i < _mesh.allvertices().size(); i++) {
		std::array<double, 3> vert;
		vert[0] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).x();
		vert[1] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).y();
		vert[2] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).z();
		verts.push_back(vert);
	}

	for (auto fp : _mesh.allfaces()) {
		auto fh = fp.first;
		auto face = fp.second;
		std::vector<uint32_t> oneface;
		for (int i = 0; i < face.getVertexHandle().size(); i++) {
			oneface.push_back(face.getVertexHandle()[i]);
		}
		faces.push_back(oneface);
	}

	int numVertices = verts.size();

	FILE* fp = fopen(_OutputFile.c_str(), "w");

	std::vector<std::string> s = SplitFileName(_OutputFile);
	s.resize(3);

	if (!fp) {
		fprintf(stdout, "Unable to open file '%s'\n", _OutputFile.c_str());
		return false;
	}

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "%s, Created by labb \n", s[1].c_str());
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	// write mesh vertices
	fprintf(fp, "POINTS %d double\n", numVertices);
	for (std::size_t i = 0; i < verts.size(); i++)
		fprintf(fp, "%.16g %.16g %.16g\n", verts[i][0], verts[i][1], verts[i][2]);
	fprintf(fp, "\n");

	int numElements = faces.size();
	//int totalNumInt = numElements * 5;
	int totalNumInt = 0;
	for (std::size_t i = 0; i < faces.size(); i++) {
		if (faces[i].size() == 3) totalNumInt += 4;
		else if (faces[i].size() == 4) totalNumInt += 5;
	}

	for (std::size_t i = 0; i < constrains.size(); ++i) {
		++numElements;
		totalNumInt += 3;
	}

	// print vertex indices in ascii or binary
	fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);
	for (std::size_t i = 0; i < constrains.size(); ++i) {
		fprintf(fp, "%d", 2);
		fprintf(fp, " %d", constrains[i][0]);
		fprintf(fp, " %d", constrains[i][1]);
		fprintf(fp, "\n");
	}

	for (std::size_t i = 0; i < faces.size(); i++) {
		fprintf(fp, "%ld", faces[i].size());
		for (int j = 0; j < faces[i].size(); j++)
			fprintf(fp, " %d", faces[i][j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// print element types in ascii or binary
	fprintf(fp, "CELL_TYPES %d\n", numElements);
	for (std::size_t i = 0; i < numElements - faces.size(); i++) {
		fprintf(fp, "%d\n", 3);
	}
	for (std::size_t i = 0; i < faces.size(); i++) {
		if (faces[i].size() == 3) fprintf(fp, "%d\n", 5);
		else if (faces[i].size() == 4) fprintf(fp, "%d\n", 9);
	}
	fclose(fp);
	return true;
}

void MeshKernel::IO::ReOrderVertexHandle(const SurfaceMesh& _mesh) {
	auto allv = _mesh.allvertices();
	int idx = 0;
	for (auto v : allv) {
		reorderedvh_.push_back(v.first);
		newvh_[v.first] = idx++;
	}
}


MeshKernel::VolumeMesh MeshKernel::IO::ReadMeshFileFromStr(const std::string& data) {
	std::stringstream ss(data);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> surface_faces;
	std::vector<std::vector<iGameVertexHandle>> cells;
	std::vector<iGameVertexHandle> cell;
	std::string str;
	enum State
	{
		USELESS = 1, VERTEX, FACE, TET, EDGE
	}state = USELESS;
	while (std::getline(ss, str)) {
		int len = str.length();
		std::vector<std::string>info;
		std::string s;
		for (int i = 0; i < len; i++) {
			if (str[i] == '#')break;
			if (str[i] == ' ') {
				if (s.length() > 0)
					info.push_back(s);
				s = "";
			}
			else {
				s.push_back(str[i]);
			}
		}
		if (s.length() > 0)
			info.push_back(s);
		if (info.size() == 0)continue;
		if (info[0] == "MeshVersionFormatted") {
			state = USELESS;
		}
		else if (info[0] == "Dimension") {
			std::getline(ss, str);
			state = USELESS;
		}
		else if (info[0] == "Vertices") {
			std::getline(ss, str);
			state = VERTEX;
		}
		else if (info[0] == "Tetrahedra") {
			std::getline(ss, str);
			state = TET;
		}
		else if (info[0] == "Triangles") {
			std::getline(ss, str);
			state = FACE;
		}
		else if (info[0] == "Edges") {
			std::getline(ss, str);
			state = EDGE;
		}
		else if (info[0] == "End") {
			state = USELESS;
		}
		else {
			if (state == USELESS) {
				continue;
			}
			else if (state == VERTEX) {
				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
					, std::stod(info[1]), std::stod(info[2])));
			}
			else if (state == TET) {
				cell = std::vector<MeshKernel::iGameVertexHandle>{
						(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
				};
				cells.push_back(cell);
				cell.clear();
			}
			else if (state == FACE) {
				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
				});
			}
		}
	}
	
	// 把模型位置移到[-1,1]范围内，中心点移到原点
    double minX = DBL_MAX; double maxX = -DBL_MAX;
    double minY = DBL_MAX; double maxY = -DBL_MAX;
    double minZ = DBL_MAX; double maxZ = -DBL_MAX;
	for (auto& v : vertices) {
		minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
		minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
		minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
	}
	double centerX = (minX + maxX) / 2;
	double centerY = (minY + maxY) / 2;
	double centerZ = (minZ + maxZ) / 2;
	double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
	for (auto& v : vertices) {
		v.x() -= centerX;
		v.y() -= centerY;
		v.z() -= centerZ;
		v = v / radius;
	}

	return VolumeMesh(vertices, cells);
}

MeshKernel::VolumeMesh MeshKernel::IO::ReadMeshFile(const std::string& _InputFile) {
//	FILE* fp = fopen(_InputFile.c_str(), "r");
//	char str[100];
//	enum State
//	{
//		USELESS = 1, VERTEX, FACE, TET, EDGE, HEX
//	}state = USELESS;
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> surface_faces;
//	std::vector<std::vector<iGameVertexHandle>> cells;
//	std::vector<iGameVertexHandle> cell;
//	while (fscanf(fp, "%[^\n]\n", str) != EOF) {
//		int barLen = strlen(str);
//		std::vector<std::string>info;
//		std::string s;
//		for (int i = 0; i < barLen; i++) {
//			if (str[i] == '#')break;
//			if (str[i] == ' ') {
//				if (s.length() > 0)
//					info.push_back(s);
//				s = "";
//			}
//			else {
//				s.push_back(str[i]);
//			}
//		}
//		if (s.length() > 0)
//			info.push_back(s);
//		if (info.size() == 0)continue;
//		if (info[0] == "MeshVersionFormatted") {
//			state = USELESS;
//		}
//		else if (info[0] == "Dimension") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = USELESS;
//		}
//		else if (info[0] == "Vertices") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = VERTEX;
//		}
//		else if (info[0] == "Tetrahedra") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = TET;
//            faces_num = 4;
//		}
//		else if (info[0] == "Triangles") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = FACE;
//		}
//		else if (info[0] == "Edges") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = EDGE;
//		}
//        else if(info[0] == "Hexahedra") {
//            fscanf(fp, "%[^\n]\n", str);
//            state = HEX;
//            faces_num = 6;
//        }
//		else if (info[0] == "End") {
//			state = USELESS;
//		}
//		else {
//			if (state == USELESS) {
//				continue;
//			}
//			else if (state == VERTEX) {
//				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
//					, std::stod(info[1]), std::stod(info[2])));
//			}
//			else if (state == TET) {
//				cell = std::vector<MeshKernel::iGameVertexHandle>{
//						(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//				};
//				cells.push_back(cell);
//				cell.clear();
//			}
//            else if (state == HEX) {
//                cell = std::vector<MeshKernel::iGameVertexHandle>{
//                        (iGameVertexHandle)(std::stoi(info[0]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[1]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[2]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[3]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[4]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[5]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[6]) - 1)
//                        , (iGameVertexHandle)(std::stoi(info[7]) - 1)
//                };
//                cells.push_back(cell);
//                cell.clear();
//            }
//			else if (state == FACE) {
//				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//			}
//		}
//	}
//    printf("vertices: %d, cells: %d\n", vertices.size(), cells.size());
//	return VolumeMesh(vertices, cells);


    string header;

    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> cells;

    ifstream stream(_InputFile, ifstream::in | ifstream::binary);

    int precision;
    int dimension;

    MeshKernel::VolumeMesh mesh;

    while (stream.good()) {
        // Read a line
        stream >> header;
        if (header.compare("MeshVersionFormatted") == 0) {
            stream >> precision;
            //HL_ASSERT_LOG(stream >> precision, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
        } else if (header.compare("Dimension") == 0) {
            stream >> dimension;
            //HL_ASSERT_LOG(stream >> dimension, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
        } else if (header.compare("Vertices") == 0) {
            int vertices_count;
            /*HL_ASSERT_LOG(stream >> vertices_count, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
            HL_LOG("[Loader] Reading %d vertices...\n", vertices_count);*/
            stream >> vertices_count;
            printf("[Read HexMesh] Reading %d vertices...\n", vertices_count);
            vertices.reserve(vertices_count);
            for (int i = 0; i < vertices_count; ++i) {
                iGameVertex v;
                float x;
                //HL_ASSERT_LOG(stream >> v.x() >> v.y() >> v.z() >> x, "ERROR: malformed mesh file. Unexpected vertex data format at vert %i.\n", i);
                stream >> v.x() >> v.y() >> v.z() >> x;
                vertices.push_back(v);
            }
        } else if (header.compare("Quadrilaterals") == 0 || header.compare("Quads") == 0) {
            int quads_count;
            stream >> quads_count;
            for (int i = 0; i < quads_count; ++i) {
                int idx[4];
                int x;
                stream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> x;
            }
        } else if (header.compare("Tetrahedra") == 0) {
            int tets_count;
            /*HL_ASSERT_LOG(stream >> hexas_count, "ERROR: malformed mesh file. Unexpected tag after hexahedras tag.\n");
            HL_LOG("[Loader] Reading %d hexas...\n", hexas_count);*/
            stream >> tets_count;
            printf("[Read HexMesh] Reading %d hexas...\n", tets_count);
            cells.reserve(tets_count);
            for (int h = 0; h < tets_count; ++h) {
                int idx[4];
                std::vector<iGameVertexHandle> vhs;
                int x;
                // irrational to rational vertex ordering!
                //stream >> idx[0] >> idx[1] >> idx[3] >> idx[2] >> idx[4] >> idx[5] >> idx[7] >> idx[6] >> x;// HexaLab
                stream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> x;// iGame
                for (int i = 0; i < 4; ++i) {
                    vhs.push_back(iGameVertexHandle(idx[i] - 1));
                }
                cells.push_back(vhs);
            }
        } else if (header.compare("Hexahedra") == 0) {
            int hexas_count;
            /*HL_ASSERT_LOG(stream >> hexas_count, "ERROR: malformed mesh file. Unexpected tag after hexahedras tag.\n");
            HL_LOG("[Loader] Reading %d hexas...\n", hexas_count);*/
            stream >> hexas_count;
            printf("[Read HexMesh] Reading %d hexas...\n", hexas_count);
            cells.reserve(hexas_count);
            for (int h = 0; h < hexas_count; ++h) {
                int idx[8];
                std::vector<iGameVertexHandle> vhs;
                int x;
                // irrational to rational vertex ordering!
                //stream >> idx[0] >> idx[1] >> idx[3] >> idx[2] >> idx[4] >> idx[5] >> idx[7] >> idx[6] >> x;// HexaLab
                stream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >> idx[6] >> idx[7] >> x;// iGame
                for (int i = 0; i < 8; ++i) {
                    vhs.push_back(iGameVertexHandle(idx[i] - 1));
                }
                cells.push_back(vhs);
            }
        } else if (header.compare("Triangles") == 0) {
            int tri_count;
            stream >> tri_count;
            for (int i = 0; i < tri_count; ++i) {
                int idx[3];
                int x;
                stream >> idx[0] >> idx[1] >> idx[2] >> x;
            }
        } else if (header.compare("Edges") == 0) {
            int edge_count;
            stream >> edge_count;
            for (int i = 0; i < edge_count; ++i) {
                int idx[2];
                int x;
                stream >> idx[0] >> idx[1] >> x;
            }
        } else if (header.compare("Corners") == 0) {
            int corner_count;
            stream >> corner_count;
            for (int i = 0; i < corner_count; ++i) {
                int c;
                stream >> c;
            }
        } else if (header.compare("End") == 0) {
            break;
        } else if (header[0] == '#') {
            stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else {
            printf("[Loader] Unexpected header \"%s\"\n ", header.c_str());
            printf("[Loader] ERROR: malformed mesh file.Unexpected header tag.\n");
            //return mesh;
        }
    }

    // Make sure at least vertex and hexa index data was read
    if (vertices.size() == 0) {
        printf("ERROR: mesh does not contain any vertex!\n");
        return mesh;
    } else if (cells.size() == 0) {
        printf("ERROR: mesh does not contain any thetra or hexa!\n");
        return mesh;
    }

    // 把模型位置移到[-1,1]范围内，中心点移到原点
    double minX = DBL_MAX; double maxX = -DBL_MAX;
    double minY = DBL_MAX; double maxY = -DBL_MAX;
    double minZ = DBL_MAX; double maxZ = -DBL_MAX;
    for (auto& v : vertices) {
        minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
        minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
        minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
    }
    double centerX = (minX + maxX) / 2;
    double centerY = (minY + maxY) / 2;
    double centerZ = (minZ + maxZ) / 2;
    double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
    for (auto& v : vertices) {
        v.x() -= centerX;
        v.y() -= centerY;
        v.z() -= centerZ;
        v = v / radius;
    }

    return VolumeMesh(vertices, cells);


}

MeshKernel::VolumeMesh MeshKernel::IO::ReadVtkFile_Volume(const std::string& _InputFile)
{
    /*
    ZQF`s Order
           v
    3----------2
    |\     ^   |\
    | \    |   | \
    |  \   |   |  \
    |   7------+---6
    |   |  +-- |-- | -> u : (+) is the barycenter of the cube
    0---+---\--1   |
     \  |    \  \  |
      \ |     \  \ |
       \|      w  \|
        4----------5
    face = { {0,3,2,1},{0,4,7,3},{0,1,5,4},{4,5,6,7},{1,2,6,5},{2,3,7,6} }
*/
    /*
        VTK File Format Order
               v
        7----------6
        |\     ^   |\
        | \    |   | \
        |  \   |   |  \
        |   4------+---5
        |   |  +-- |-- | -> u : (+) is the barycenter of the cube
        3---+---\--2   |
         \  |    \  \  |
          \ |     \  \ |
           \|      w  \|
            0----------1
    */

    std::ifstream hexfile(_InputFile.c_str());

    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> cells;

    if (hexfile.is_open())
    {
        std::cout << "Reading " << _InputFile << " File" << std::endl;

        std::string str;
        try
        {
            do
            {
                hexfile >> str;
            } while (str != "DATASET");
            hexfile >> str;
            if (str != "UNSTRUCTURED_GRID")
            {
                //should stop reading, but not implement here
                1 + 1 == 2;
            }
            do
            {
                hexfile >> str;
            } while (str != "POINTS");

            size_t nv, ne, num_ne_type;

            hexfile >> nv >> str;

            vertices.resize(nv);
            for (size_t i = 0; i < nv; i++)
            {
                hexfile >> vertices[i].x() >> vertices[i].y() >> vertices[i].z();
            }

            hexfile >> str >> ne >> num_ne_type;

            cells.reserve(ne);
            std::vector<iGameVertexHandle> vhs_per_cell;
            vhs_per_cell.resize(8);
            int index[8] = { 4, 5, 1, 0, 7, 6, 2, 3 };
            int tmp = 0;
            for (size_t i = 0; i < ne; i++)
            {
                //处理掉第一个numPoints Flag
                hexfile >> tmp; // 8

                for (int j = 0; j < 8; j++)
                {
                    hexfile >> tmp;
                    vhs_per_cell[index[j]] = iGameVertexHandle(tmp);
                }
                cells.push_back(vhs_per_cell);
            }
        }
        catch (...)
        {
            hexfile.close();
        }
        hexfile.close();
        std::cout << "Reading " << " Finished" << std::endl;
    }

    // 把模型位置移到[-1,1]范围内，中心点移到原点
    double minX = DBL_MAX; double maxX = -DBL_MAX;
    double minY = DBL_MAX; double maxY = -DBL_MAX;
    double minZ = DBL_MAX; double maxZ = -DBL_MAX;
    for (auto& v : vertices) {
        minX = std::min(minX, v.x()), maxX = std::max(maxX, v.x());
        minY = std::min(minY, v.y()), maxY = std::max(maxY, v.y());
        minZ = std::min(minZ, v.z()), maxZ = std::max(maxZ, v.z());
    }
    double centerX = (minX + maxX) / 2;
    double centerY = (minY + maxY) / 2;
    double centerZ = (minZ + maxZ) / 2;
    double radius = std::max(maxX - minX, std::max(maxY - minY, maxZ - minZ)) / 2;
    for (auto& v : vertices) {
        v.x() -= centerX;
        v.y() -= centerY;
        v.z() -= centerZ;
        v = v / radius;
    }

    auto mesh = MeshKernel::VolumeMesh(vertices, cells);

    return mesh;
}

// 保存体结构，目前只支持保存四面体
bool MeshKernel::IO::WriteMeshFile(const VolumeMesh& _mesh, const std::string& _OutputFile) {

    std::ofstream off(_OutputFile.c_str(), std::ios::out);

    if (!off.good()) {
        std::cerr << "Error: Could not open file " << _OutputFile << " for writing!" << std::endl;
        off.close();
        return false;
    }

    int nbv, nbt;
    nbv = _mesh.vsize();
    nbt = _mesh.CellSize();
    vector<double> vertices(3 * nbv);
    vector<int> tet(4 * nbt);
    for (auto& vp : _mesh.allvertices()) {
        auto& vh = vp.first;
        auto& v = vp.second;
        int idx = vh.idx();
        vertices[idx * 3] = v.x();
        vertices[idx * 3 + 1] = v.y();
        vertices[idx * 3 + 2] = v.z();
    }
    for (auto& cp : _mesh.allcells()) {
        auto& ch = cp.first;
        auto& c = cp.second;
        int idx = ch.idx();
        vector<VH> vex = c.getVertexHandle();
        tet[idx * 4] = vex[0];
        tet[idx * 4 + 1] = vex[1];
        tet[idx * 4 + 2] = vex[2];
        tet[idx * 4 + 3] = vex[3];
    }

    // Write header
    off << "MeshVersionFormatted 1" << std::endl;
    off << "Dimension 3" << std::endl;
    off << "Vertices" << std::endl;
    off << nbv << std::endl;

    for (int i = 0; i < nbv; i++) {
        off << vertices[3 * i] << " " << vertices[3 * i + 1] << " " << vertices[3 * i + 2] << " " << "-1" << std::endl;
    }
    off << "Tetrahedra" << std::endl;
    off << nbt << std::endl;
    for (int i = 0; i < nbt; i++) {
        off << tet[4 * i] + 1 << " " << tet[4 * i + 1] + 1 << " " << tet[4 * i + 2] + 1<< " " << tet[4 * i + 3] + 1<< " " << "1" << std::endl;
    }
    off << "End" << std::endl;
    off.close();
	return true;
}

