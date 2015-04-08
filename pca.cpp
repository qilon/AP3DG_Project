#include "pca.h"
//=============================================================================
PCA::PCA()
{

}
//=============================================================================
PCA::PCA(string _pca_filename_url)
{
	read(_pca_filename_url);
	initAlphas();
}
//=============================================================================
PCA::PCA(int _n_meshes, string _ply_models_url_preffix)
{
	degrees_freedom = _n_meshes - 1;

	initAlphas();

	// Starting counting time
	time_t tstart, tend;
	tstart = time(0);

	// Reading meshes:
	MyMesh *meshes = new MyMesh[_n_meshes];

	for (int iMesh = 0; iMesh < _n_meshes; iMesh++){
		string index;
		stringstream convert;
		convert << iMesh;
		index = convert.str();
		cout << "Reading mesh #" + index << endl;
		if (!OpenMesh::IO::read_mesh(meshes[iMesh], _ply_models_url_preffix + index + ".ply"))
		{
			std::cerr << "Cannot read mesh #" + index << std::endl;
		}
	}

	// Building matrix S:
	cout << "Building matrix S..." << endl;
	const int nVert = meshes[0].n_vertices();
	Eigen::MatrixXf S(3 * nVert, _n_meshes);
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++) {
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
	Eigen::MatrixXf centS(3 * nVert, _n_meshes);
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++) {
		centS.col(iMesh) = S.col(iMesh) - meanS;
	}

	// Building covariance matrix C:
	cout << "Building covariance matrix C..." << endl;
	Eigen::MatrixXf C = centS * centS.transpose();

	// Computing eigenvectors:
	cout << "Computing eigenvectors..." << endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(C);

	eigen_vectors = solver.eigenvectors().real();
	eigen_values = solver.eigenvalues().real();

	// Deleting meshes
	delete meshes;
	meshes = nullptr;

	// Finishing counting time
	tend = time(0);

	cout << "Elapsed time : " << difftime(tend, tstart) / 60 << " minute(s)." << endl;
}
//=============================================================================
PCA::~PCA()
{
}
//=============================================================================
void PCA::read(string _pca_filename_url)
{
	MatrixXf::Index evectors_rows = 0;
	MatrixXf::Index evectors_cols = 0;
	MatrixXf::Index evalues_rows = 0;
	MatrixXf::Index evalues_cols = 0;

	ifstream in(_pca_filename_url, ios::in | std::ios::binary);

	in.read((char*)(&degrees_freedom), sizeof(int));

	in.read((char*)(&evectors_rows), sizeof(MatrixXf::Index));
	in.read((char*)(&evectors_cols), sizeof(MatrixXf::Index));
	in.read((char*)(&evalues_rows), sizeof(MatrixXf::Index));
	in.read((char*)(&evalues_cols), sizeof(MatrixXf::Index));

	eigen_vectors.resize(evectors_rows, evectors_cols);
	eigen_values.resize(evalues_rows, evalues_cols);

	in.read((char *)eigen_vectors.data(), 
		evectors_rows*evectors_cols*sizeof(MatrixXf::Scalar));
	in.read((char *)eigen_values.data(),
		evalues_rows*evalues_cols*sizeof(MatrixXf::Scalar));

	in.close();
}
//=============================================================================
void PCA::write(string _pca_filename_url)
{
	MatrixXf::Index evectors_rows = eigen_vectors.rows();
	MatrixXf::Index evectors_cols = eigen_vectors.cols();
	MatrixXf::Index evalues_rows = eigen_values.rows();
	MatrixXf::Index evalues_cols = eigen_values.cols();

	ofstream out(_pca_filename_url, ios::out | ios::binary | ios::trunc);

	out.write((char*)(&degrees_freedom), sizeof(int));

	out.write((char*)(&evectors_rows), sizeof(MatrixXf::Index));
	out.write((char*)(&evectors_cols), sizeof(MatrixXf::Index));
	out.write((char*)(&evalues_rows), sizeof(MatrixXf::Index));
	out.write((char*)(&evalues_cols), sizeof(MatrixXf::Index));

	out.write((char*)eigen_vectors.data(), 
		evectors_rows*evectors_cols*sizeof(MatrixXf::Scalar));
	out.write((char*)eigen_values.data(),
		evalues_rows*evalues_cols*sizeof(MatrixXf::Scalar));

	out.close();
}
//=============================================================================
void PCA::initAlphas()
{
	alphas = VectorXf::Zero(degrees_freedom);
	alphas(0) = 1;
}
//=============================================================================