#include "GLViewer.h"
//=============================================================================
// PCA PARAMETERS
const static int NUM_MESHES = 71;
const static string PLY_MODELS_URL_PREFFIX = "_models/scapecomp/mesh";
const static string PCA_RESULT_URL = "./_data/new_pca_result.dat";
const static string FEATURES_URL = "./_data/features.dat";

// VIEWER TEST PARAMETERS
const static char* PLY_FILENAME = "./_models/scapecomp/mesh0.ply";

//=============================================================================
int main(int argc, char **argv)
{
	
	/* Generate the eigen vectors and eigen values of the covariance matrix for 
	the 71 meshes and keeps the 70 highest ones */
	//PCA pca = PCA(NUM_MESHES, PLY_MODELS_URL_PREFFIX);
	//pca.write(PCA_RESULT_URL);
	//pca.writeFeatures(NUM_MESHES, PLY_MODELS_URL_PREFFIX, FEATURES_URL);

	// PCA pca = PCA(PCA_RESULT_URL, FEATURES_URL);
	// pca.writeFeatures(NUM_MESHES, PLY_MODELS_URL_PREFFIX, FEATURES_URL);
	
	// Mesh viewer
	GLViewer viewer;
	viewer.initialize(&argc, argv);

	// Reads a mesh to get the vertices connections
	MyMesh mesh;
	loadMesh(mesh, PLY_FILENAME);
	viewer.setMesh(mesh);

	// Loads PCA info
	viewer.loadPCA(PCA_RESULT_URL, FEATURES_URL);

	// Loop viewer
	viewer.run();

	return 0;
}
//=============================================================================