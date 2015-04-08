#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "convert.h"
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

private:
	int degrees_freedom;
	MatrixXf eigen_vectors;
	MatrixXf eigen_values;
	VectorXf alphas;

	void initAlphas();
};
//=============================================================================