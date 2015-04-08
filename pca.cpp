#include "pca.h"
//=============================================================================

void pca(int nMeshes, string ply_models_url_preffix, string pca_result_url)
{
	// Starting counting time
	time_t tstart, tend;
	tstart = time(0);
	
	// Reading meshes:
	MyMesh *meshes = new MyMesh[nMeshes];

	for (int iMesh = 0; iMesh < nMeshes; iMesh++){
		string index;
		stringstream convert;
		convert << iMesh;
		index = convert.str();
		cout << "Reading mesh #" + index << endl;
		if (!OpenMesh::IO::read_mesh(meshes[iMesh], ply_models_url_preffix + index + ".ply"))
		{
			std::cerr << "Cannot read mesh #" + index << std::endl;
		}
	}

	// Building matrix S:
	cout << "Building matrix S..." << endl;
	const int nVert = meshes[0].n_vertices();
	Eigen::MatrixXf S(3 * nVert, nMeshes);
	for (int iMesh = 0; iMesh < nMeshes; iMesh++) {
		Eigen::MatrixXf cloud = mesh2EigenMatrix(meshes[iMesh]);
		for (int iVert = 0; iVert < nVert; iVert++) {
			S(3 * iVert, iMesh) = cloud(0, iVert);
			S(3 * iVert + 1, iMesh) = cloud(1, iVert);
			S(3 * iVert + 2, iMesh) = cloud(2, iVert);
		}
	}
	
	// Computing mean:
	cout << "Building mean..." << endl;
	Eigen::VectorXf meanS = S.rowwise().mean();
	
	// Computing centered S:
	cout << "Computing centered S..." << endl;
	Eigen::MatrixXf centS(3 * nVert, nMeshes);
	for (int iMesh = 0; iMesh < nMeshes; iMesh++) {
		centS.col(iMesh) = S.col(iMesh) - meanS;
	}

	// Building covariance matrix C:
	cout << "Building covariance matrix C..." << endl;
	Eigen::MatrixXf C = centS * centS.transpose();
	
	// Computing eigenvectors:
	cout << "Computing eigenvectors..." << endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(C);
	Eigen::MatrixXf eVectors = solver.eigenvectors().real();
	Eigen::VectorXf eValues = solver.eigenvalues().real();

	// Deleting meshes
	delete meshes;
	meshes = nullptr;

	// Finishing counting time
	tend = time(0);
	
	// Writing result in file
	ofstream file(pca_result_url);
	if (file.is_open())
	{
		file << "Elapsed time : " << difftime(tend, tstart)/60 << " minute(s)." << endl;
		file << "Eigen vectors: " << endl << eVectors << endl;
		file << "Eigen values: " << endl << eValues << endl;
	}
}
//=============================================================================