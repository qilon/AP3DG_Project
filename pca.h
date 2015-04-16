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
#define RIGHT_ELBOW 8616
#define LEFT_ELBOW 8793
#define RIGHT_KNUCKLE 5165
#define LEFT_KNUCKLE 5736
#define RIGHT_SHOULDER 10721
#define LEFT_SHOULDER 10992
#define RIGHT_AXILLA 9784
#define LEFT_AXILLA 9919
#define RIGHT_INSIDE_ELBOW 8344
#define LEFT_INSIDE_ELBOW 8500
#define LEFT_WAIST 7796
#define RIGHT_WAIST 7665
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
};
//=============================================================================