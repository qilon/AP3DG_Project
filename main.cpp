#include "GLViewer.h"
//=============================================================================
#define COMPUTE_PCA 0				// Compute PCA model 
#define COMPUTE_FEATURES 1			// Compute features mapping 
#define RUN_VIEWER 2				// Run viewer
#define GENERATE_NOISE 3			// Generate noisy data

#define BODY_DATASET		0
#define NOISY_BODY_DATASET	1
#define SCAPE_DATASET		2
//=============================================================================
/*
 * CHANGE THIS VALUES TO RUN DIFFERENT TASKS WITH DIFFERENT DATASETS
 */
const static int DATASET = BODY_DATASET;
const static int TASK = RUN_VIEWER;

//=============================================================================
int main(int argc, char **argv)
{
	PCA pca;
	GLViewer viewer;
	string pca_model_url, features_url;
	string models_url_preffix, models_url_suffix;
	string noise_result_url_preffix;
	char* mesh_filename;
	int num_meshes, mesh_first_idx;

	const static string BODY_FEATURES_URL = "./_data/body_data_features_no_zeros.txt";

	switch (DATASET)
	{
	case BODY_DATASET:
		pca_model_url = "./_data/body_pca_model_no_zeros.dat";
		features_url = "./_data/body_features_no_zeros.dat";
		models_url_preffix = "_models/body_data_no_zeros/s";
		models_url_suffix = "p0.obj";
		num_meshes = 108;
		mesh_first_idx = 1;
		mesh_filename = "./_models/body_data_no_zeros/s58p0.obj";
		noise_result_url_preffix = "_models/body_data_no_zeros_noise/s";
		break;

	case NOISY_BODY_DATASET:
		pca_model_url = "./_data/body_pca_model_no_zeros_noise.dat";
		features_url = "./_data/body_features_no_zeros_noise.dat";
		models_url_preffix = "_models/body_data_no_zeros_noise/s";
		models_url_suffix = "p0.obj";
		num_meshes = 108;
		mesh_first_idx = 1;
		mesh_filename = "./_models/body_data_no_zeros/s58p0.obj";
		noise_result_url_preffix = "_models/body_data_no_zeros_noise/s";
		break;

	case SCAPE_DATASET:
		pca_model_url = "./_data/new_pca_result.dat";
		features_url = "./_data/new_features.dat";
		models_url_preffix = "_models/scapecomp/mesh";
		models_url_suffix = ".ply";
		num_meshes = 71;
		mesh_first_idx = 0;
		mesh_filename = "./_models/scapecomp/mesh0.ply";
		noise_result_url_preffix = "_models/scapecomp_noise/mesh";
		break;
	default:
		break;
	}

	switch (TASK)
	{
	case COMPUTE_PCA:	// Compute PCA model 
		pca = PCA(num_meshes, models_url_preffix, models_url_suffix,
			mesh_first_idx, pca_model_url);
		break;

	case COMPUTE_FEATURES:	// Compute features mapping
		pca = PCA(pca_model_url);
		if (DATASET==SCAPE_DATASET)
		{
			pca.computeFeatures(num_meshes, models_url_preffix, features_url);
		}
		else
		{
			pca.computeFeatures(BODY_FEATURES_URL, models_url_preffix,
				models_url_suffix, mesh_first_idx, features_url);
		}
		break;

	case RUN_VIEWER:	// Run mesh viewer
		viewer.initialize(&argc, argv);
		viewer.loadMesh(mesh_filename);	// Reads a mesh to get the vertices connections
		viewer.loadPCA(pca_model_url, features_url); // Loads PCA info
		viewer.run();
		viewer.destroy();
		break;

	case GENERATE_NOISE:	// Add noise to meshes
		addNoise(models_url_preffix, models_url_suffix, mesh_first_idx,
			num_meshes, noise_result_url_preffix);
		break;

	default:
		break;
	}

	return 0;
}
//=============================================================================