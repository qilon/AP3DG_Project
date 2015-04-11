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
using namespace std;
using namespace Eigen;
//=============================================================================
class PCA
{
public:
	PCA();
	PCA(int _n_meshes, string _ply_models_url_preffix);
	PCA(string _pca_filename_url);
	~PCA();
	void read(string _pca_filename_url);
	void write(string _pca_filename_url);
	void updateMesh(MyMesh& _mesh);
	void editAlpha(int idx, float new_value);
	float getAlpha(int idx);

private:
	int degrees_freedom;
	int vector_size;
	MatrixXf eigen_vectors;
	VectorXf eigen_values;
	VectorXf mean_model;
	VectorXf first_model;
	VectorXf alphas;

	void initAlphas();
	void computePCA(MyMesh* meshes, int _n_meshes);
};
//=============================================================================