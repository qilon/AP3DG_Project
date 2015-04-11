#include "pca.h"
//=============================================================================
PCA::PCA()
{
	degrees_freedom = 0;
	vector_size = 0;
}
//=============================================================================
PCA::PCA(string _pca_filename_url)
{
	read(_pca_filename_url);
}
//=============================================================================
PCA::PCA(int _n_meshes, string _ply_models_url_preffix)
{
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

	// Computing PCA
	computePCA(meshes, _n_meshes);

	// Initialize coefficients
	initAlphas();

	// Deleting meshes
	delete[] meshes;
	meshes = nullptr;
}
//=============================================================================
PCA::~PCA()
{
}
//=============================================================================
void PCA::read(string _pca_filename_url)
{
	ifstream in(_pca_filename_url, ios::in | std::ios::binary);

	in.read((char*)(&degrees_freedom), sizeof(int));
	in.read((char*)(&vector_size), sizeof(int));

	cout << "Degrees of freedom: " << degrees_freedom << endl;
	cout << "Eigenvectors rows: " << vector_size << endl; 

	eigen_vectors.resize(vector_size, degrees_freedom);
	eigen_values.resize(degrees_freedom);
	mean_model.resize(vector_size);

	in.read((char *)eigen_vectors.data(),
		vector_size*degrees_freedom*sizeof(MatrixXf::Scalar));
	in.read((char *)eigen_values.data(),
		degrees_freedom*sizeof(VectorXf::Scalar));
	in.read((char*)mean_model.data(),
		vector_size*sizeof(VectorXf::Scalar));

	in.close();

	initAlphas();

	// First model:
	MyMesh first;
	const char* filename = "./_models/scapecomp/mesh0.ply";
	loadMesh(first, filename, false);
	Eigen::MatrixXf cloud = mesh2EigenMatrix(first);
	int nVert = cloud.cols();
	cout << vector_size << endl;
	first_model.resize(vector_size);
	for (int iVert = 0; iVert < nVert; iVert++) {
		first_model(3 * iVert) = cloud(0, iVert);
		first_model(3 * iVert + 1) = cloud(1, iVert);
		first_model(3 * iVert + 2) = cloud(2, iVert);
	}
}
//=============================================================================
void PCA::write(string _pca_filename_url)
{
	ofstream out(_pca_filename_url, ios::out | ios::binary | ios::trunc);

	out.write((char*)(&degrees_freedom), sizeof(int));
	out.write((char*)(&vector_size), sizeof(int));

	cout << "Degrees of freedom: " << degrees_freedom << endl;
	cout << "Eigenvectors rows: " << vector_size << endl;

	out.write((char*)eigen_vectors.data(),
		vector_size*degrees_freedom*sizeof(MatrixXf::Scalar));
	out.write((char*)eigen_values.data(),
		degrees_freedom*sizeof(VectorXf::Scalar));
	out.write((char*)mean_model.data(),
		vector_size*sizeof(VectorXf::Scalar));

	out.close();
}
//=============================================================================
void PCA::computePCA(MyMesh* meshes, int _n_meshes)
{
	int nVert = meshes[0].n_vertices();

	// Building matrix S:
	cout << "Building matrix S..." << endl;
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
	VectorXf meanS = S.rowwise().mean();

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

	MatrixXf new_eigen_vectors = solver.eigenvectors().real();
	VectorXf new_eigen_values = solver.eigenvalues().real();

	// Setting mean model, degrees of freedom and size vector
	mean_model = meanS;
	degrees_freedom = _n_meshes - 1;
	vector_size = new_eigen_vectors.rows();
	eigen_vectors.resize(vector_size, degrees_freedom);
	eigen_values.resize(degrees_freedom);

	// Taking the highest eigenvalues and their corresponding eigenvectors 
	float min_eigen_value = new_eigen_values.cwiseAbs().minCoeff() - 1.f;
	VectorXf::Index max_idx = 0;
	for (int i = 0; i < degrees_freedom; i++)
	{
		new_eigen_values.cwiseAbs().maxCoeff(&max_idx);
		eigen_vectors.col(i) = new_eigen_vectors.col(max_idx);
		eigen_values(i) = new_eigen_values(max_idx);
		new_eigen_values(max_idx) = min_eigen_value;
	}
}
//=============================================================================
void PCA::updateMesh(MyMesh& _mesh)
{
	VectorXf new_model = first_model;

	// cout << alphas.size() << endl;
	// cout << eigen_vectors.rows() << endl;
	// cout << eigen_vectors.cols() << endl;

	for (int i = 0; i < degrees_freedom; i++)
	{
		new_model += alphas(i) * eigen_vectors.col(i);
	}

	MyMesh::ConstVertexIter vIt1(_mesh.vertices_begin());
	MyMesh::ConstVertexIter vIt2(_mesh.vertices_end());

	for (size_t i = 0; vIt1 != vIt2; vIt1++, i++) {
		float x = new_model(3 * i);
		float y = new_model(3 * i + 1);
		float z = new_model(3 * i + 2);
		_mesh.point(*vIt1)[0] = x;
		_mesh.point(*vIt1)[1] = y;
		_mesh.point(*vIt1)[2] = z;
	}
}
//=============================================================================
void PCA::initAlphas()
{
	alphas = VectorXf::Zero(degrees_freedom);
	// alphas = eigen_values;
}
//=============================================================================
void PCA::editAlpha(int idx, float new_value)
{
	alphas(idx) = new_value;
}
//=============================================================================
float PCA::getAlpha(int idx)
{
	return alphas(idx);
}
//=============================================================================