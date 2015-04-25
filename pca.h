#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "MyMesh.h"
#include "FeatureConfig.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
//=============================================================================
#define RIGHT_ELBOW 8616
#define LEFT_ELBOW 8793
#define RIGHT_BACK_WRIST 6475
#define RIGHT_FRONT_WRIST 6628
#define LEFT_BACK_WRIST 6655
#define LEFT_FRONT_WRIST 6870
#define RIGHT_KNUCKLE 5165
#define LEFT_KNUCKLE 5736
#define RIGHT_SHOULDER 10721
#define LEFT_SHOULDER 10992
#define RIGHT_AXILLA 9784
#define LEFT_AXILLA 9919
#define RIGHT_INSIDE_ELBOW 8344
#define LEFT_INSIDE_ELBOW 8500
#define RIGHT_WAIST 7665
#define LEFT_WAIST 7796
#define RIGHT_KNEE 2850
#define LEFT_KNEE 2863
#define RIGHT_INSIDE_ANKLE 1127
#define LEFT_INSIDE_ANKLE 1055
#define RIGHT_OUTSIDE_ANKLE 1082
#define LEFT_OUTSIDE_ANKLE 1101
#define RIGHT_INSIDE_KNEE 2703
#define LEFT_INSIDE_KNEE 2667
#define RIGHT_HIP 6304
#define LEFT_HIP 6017
#define PERINEUM 12492
//=============================================================================
using namespace std;
using namespace Eigen;
//=============================================================================
class PCA
{
public:
	PCA();
	PCA(int _n_meshes, string _ply_models_url_preffix, 
		string _ply_models_url_suffix, int first_index, string _pca_filename_url);
	PCA(string _pca_filename_url, string _features_filename_url);
	PCA(string _pca_filename_url);
	virtual ~PCA();
	void readPCA(string _pca_filename_url);
	void computeFeatures(string _features_data_filename_url, string _ply_models_url_preffix,
		string _ply_models_url_suffix, int first_index, string _feature_filename_url);
	void computeFeatures(int _n_meshes, string _ply_models_url_preffix, string _feature_filename_url);
	void readFeatures(string _features_filename_url);
	void updateMesh(MyMesh& _mesh);
	void PCA::editAlpha(int _idxAlpha, float value);
	VectorXf getAlphas();
	void editFeature(int idxFeature, float new_value);
	VectorXf getFeatures();
	FeatureConfig* getInitialFeatures();
	int getNFeatures();
	void reconstructMesh(MyMesh& _mesh, const VectorXi& _points_state,
		const MyMesh::Color& color);

private:
	int degrees_freedom;
	int vector_size;
	int n_features;
	MatrixXf eigen_vectors;
	VectorXf eigen_values;
	VectorXf mean_model;
	VectorXf first_model;
	VectorXf alphas;
	VectorXf features;
	MatrixXf M_feature2Alpha;
	FeatureConfig* initialFeatures;

	void computePCA(MyMesh* meshes, int _n_meshes, string _pca_filename_url);
	void writePCA(string _pca_filename_url);
	MatrixXf readFeaturesData(string _features_data_filename_url);
	void writeFeatures(string _feature_filename_url);
	void centerModel();
	void initAlphas();
	void initFeatures();
	void updateMesh(MyMesh& _mesh, const VectorXf& _v_mesh,
		const VectorXi& _points_state, const MyMesh::Color& _color);
	void updateMeshAlphas(MyMesh& _mesh, const VectorXf& _alphas,
		const VectorXi& _points_state, const MyMesh::Color& _color);
	MatrixXf pinv(const MatrixXf& _m);
};
//=============================================================================