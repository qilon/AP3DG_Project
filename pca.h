#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "MyMesh.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
//=============================================================================
#define PCA_FILE 0
#define FULL_EIGEN_FILE 1
//=============================================================================
using namespace std;
using namespace Eigen;
//=============================================================================
class PCA
{
public:
	PCA();
	PCA(int _n_meshes, string _ply_models_url_preffix);
	PCA(string _pca_filename_url, int _file_type = PCA_FILE);
	~PCA();
	void readFullEigen(string _pca_filename_url);
	void writeFullEigen(string _pca_filename_url);
	void read(string _pca_filename_url);
	void write(string _pca_filename_url);
	void updateMesh(MyMesh& _mesh);

private:
	int degrees_freedom;
	int vector_size;
	MatrixXf eigen_vectors;
	VectorXf eigen_values;
	VectorXf mean_model;
	VectorXf alphas;

	void initAlphas();
	void computePCA(MyMesh* meshes, int _n_meshes);
};
//=============================================================================