#include "GLViewer.h"
//=============================================================================
// PCA PARAMETERS
const static int NUM_MESHES = 108;
const static string MODELS_URL_PREFFIX = "_models/body_data_no_zeros/s";
const static string MODELS_URL_SUFFIX = "p0.obj";
const static string MODELS_NOISE_URL_PREFFIX = "_models/body_data_no_zeros_noise/s";
const static int MODELS_URL_FIRST_IDX = 1;
const static string PCA_MODEL_URL = "./_data/body_pca_model_no_zeros.dat";
const static string PCA_MODEL_NOISE_URL = "./_data/body_pca_model_no_zeros_noise.dat";
const static string FEATURES_URL = "./_data/body_features_no_zeros.dat";
const static string BODY_FEATURES_URL = "./_data/body_data_features_no_zeros.txt";

// VIEWER TEST PARAMETERS
const static char* MESH_FILENAME = "./_models/body_data_no_zeros/s1p0.obj";

//=============================================================================
int main(int argc, char **argv)
{
	// Add noise to meshes
	//addNoise(MODELS_URL_PREFFIX, MODELS_URL_SUFFIX, MODELS_URL_FIRST_IDX,
	//	NUM_MESHES, MODELS_NOISE_URL_PREFFIX);

	// Compute PCA model 
	//PCA pca = PCA(NUM_MESHES, MODELS_NOISE_URL_PREFFIX, MODELS_URL_SUFFIX, 
	//	MODELS_URL_FIRST_IDX, PCA_MODEL_NOISE_URL);

	// Compute features mapping
	// Load PCA info only:

	// PCA pca = PCA(PCA_MODEL_NOISE_URL);
	// Body Data:
	// pca.computeFeatures(BODY_FEATURES_URL, MODELS_URL_PREFFIX, MODELS_URL_SUFFIX, MODELS_URL_FIRST_IDX, FEATURES_URL);
	// pca.computeFeatures(BODY_FEATURES_URL, MODELS_NOISE_URL_PREFFIX, MODELS_URL_SUFFIX, MODELS_URL_FIRST_IDX, FEATURES_URL);

	// Scapecomp:
	// pca.computeFeatures(NUM_MESHES, MODELS_URL_PREFFIX, FEATURES_URL);

	// Run mesh viewer
	GLViewer viewer;
	viewer.initialize(&argc, argv);
	viewer.loadMesh(MESH_FILENAME);	// Reads a mesh to get the vertices connections
	viewer.loadPCA(PCA_MODEL_URL, FEATURES_URL); // Loads PCA info
	// viewer.loadPCA(PCA_MODEL_NOISE_URL, FEATURES_URL); // Loads PCA info for noisy meshes
	viewer.run();
	viewer.destroy();

	return 0;
}
//=============================================================================