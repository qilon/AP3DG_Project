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
	PCA(string _pca_filename_url, string _features_filename_url);
	~PCA();
	void read(string _pca_filename_url, string _features_filename_url);
	void write(string _pca_filename_url);
	void updateMesh(MyMesh& _mesh);
	void editFeature(int idxFeature, float new_value);
	float PCA::getFeature(int idxFeature);
	int PCA::getControllers();
	void PCA::writeFeatures(int _n_meshes, string _ply_models_url_preffix,
		string _feature_filename_url);

private:
	int degrees_freedom;
	int vector_size;
	int n_controllers;
	MatrixXf eigen_vectors;
	VectorXf eigen_values;
	VectorXf mean_model;
	VectorXf first_model;
	VectorXf alphas;
	VectorXf features;
	MatrixXf M_feature2Alpha;

	void initAlphas();
	void computePCA(MyMesh* meshes, int _n_meshes);
	void centerModel();
};
//=============================================================================