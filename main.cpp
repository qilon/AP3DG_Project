#include "GLViewer.h"
//=============================================================================
// PCA PARAMETERS
const static int NUM_MESHES = 111;
const static string MODELS_URL_PREFFIX = "_models/body_data/s";
const static string MODELS_URL_SUFFIX = "p0.obj";
const static int MODELS_URL_FIRST_IDX = 1;
const static string PCA_MODEL_URL = "./_data/body_pca_model.dat";
const static string FEATURES_URL = "./_data/features.dat";

// VIEWER TEST PARAMETERS
const static char* MESH_FILENAME = "./_models/body_data/s1p0.obj";

//=============================================================================
int main(int argc, char **argv)
{
	
	// Generate the eigen vectors and eigen values of the covariance matrix for 
	// the 71 meshes and keeps the 70 highest ones
	//PCA pca = PCA(NUM_MESHES, MODELS_URL_PREFFIX, MODELS_URL_SUFFIX, 
		//MODELS_URL_FIRST_IDX);
	//pca.writePCA(PCA_MODEL_URL);
	//pca.writeFeatures(NUM_MESHES, PLY_MODELS_URL_PREFFIX, FEATURES_URL);

	// PCA pca = PCA(PCA_MODEL_URL);
	// pca.computeFeatures(NUM_MESHES, PLY_MODELS_URL_PREFFIX);
	// pca.writeFeatures(FEATURES_URL);

	//PCA pca = PCA(PCA_MODEL_URL, FEATURES_URL);

	// Mesh viewer
	GLViewer viewer;
	viewer.initialize(&argc, argv);

	// Reads a mesh to get the vertices connections
	MyMesh mesh;
	loadMesh(mesh, MESH_FILENAME);
	viewer.setMesh(mesh);

	// Loads PCA info
	viewer.loadPCA(PCA_MODEL_URL, FEATURES_URL);

	// Loop viewer
	viewer.run();

	return 0;
}
//=============================================================================