#include "pca.h"
//=============================================================================
PCA::PCA()
{
	degrees_freedom = 0;
	vector_size = 0;
}
//=============================================================================
PCA::PCA(string _pca_filename_url, string _features_filename_url)
{
	read(_pca_filename_url, _features_filename_url);
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
void PCA::read(string _pca_filename_url, string _features_filename_url)
{
	// Load PCA info
	ifstream inPCA(_pca_filename_url, ios::in | std::ios::binary);

	inPCA.read((char*)(&degrees_freedom), sizeof(int));
	inPCA.read((char*)(&vector_size), sizeof(int));

	cout << "Degrees of freedom: " << degrees_freedom << endl;
	cout << "Eigenvectors rows: " << vector_size << endl; 

	eigen_vectors.resize(vector_size, degrees_freedom);
	eigen_values.resize(degrees_freedom);
	mean_model.resize(vector_size);
	first_model.resize(vector_size);

	inPCA.read((char *)eigen_vectors.data(),
		vector_size*degrees_freedom*sizeof(MatrixXf::Scalar));
	inPCA.read((char *)eigen_values.data(),
		degrees_freedom*sizeof(VectorXf::Scalar));
	inPCA.read((char*)mean_model.data(),
		vector_size*sizeof(VectorXf::Scalar));
	inPCA.read((char*)first_model.data(),
		vector_size*sizeof(VectorXf::Scalar));


	inPCA.close();

	
	// Load features
	// ifstream inFeatures(_features_filename_url, ios::in | std::ios::binary);
	ifstream inFeatures(_features_filename_url, ios::in);

	inFeatures.read((char*)(&n_controllers), sizeof(int));

	cout << "Number of controllers: " << n_controllers << endl;
	M_feature2Alpha.resize(degrees_freedom, n_controllers + 1);
	features.resize(n_controllers + 1);

	inFeatures.read((char *)M_feature2Alpha.data(),
		degrees_freedom*(n_controllers + 1)*sizeof(MatrixXf::Scalar));
	inFeatures.read((char *)features.data(),
		(n_controllers + 1)*sizeof(MatrixXf::Scalar));
	inFeatures.close();
	


	initAlphas();
	centerModel();
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
	out.write((char*)first_model.data(),
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
	first_model = S.col(0);

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
	float min_eigen_value = new_eigen_values.cwiseAbs().minCoeff() / 2.f;
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
	VectorXf new_model = mean_model;

	alphas = M_feature2Alpha * features;
	cout << alphas << endl;

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
void PCA::editFeature(int idxFeature, float new_value)
{
	features(idxFeature) = new_value;
}
//=============================================================================
float PCA::getFeature(int idxFeature)
{
	return features(idxFeature);
}
//=============================================================================
int PCA::getControllers()
{
	return n_controllers;
}
//=============================================================================
void PCA::writeFeatures(int _n_meshes, string _ply_models_url_preffix,
	string _feature_filename_url) {

	int n_controllers = 1;

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

	int nVert = mean_model.size()/3;
	cout << nVert << endl;
	MatrixXf features(MatrixXf::Ones(n_controllers + 1, _n_meshes));
	MatrixXf centS(3 * nVert, _n_meshes);
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++){
		MatrixXf model = mesh2EigenMatrix(meshes[iMesh]);
		Vector3f u = model.col(8344) - model.col(9784);
		Vector3f v = model.col(7665) - model.col(9784);
		float angle = acos(u.dot(v) / (u.norm() * v.norm()));
		features(0, iMesh) = angle;

		for (int iVert = 0; iVert < nVert; iVert++) {
			centS(3 * iVert, iMesh) = model(0, iVert) - mean_model(3 * iVert);
			centS(3 * iVert + 1, iMesh) = model(1, iVert) - mean_model(3 * iVert + 1);
			centS(3 * iVert + 2, iMesh) = model(2, iVert) - mean_model(3 * iVert + 2);
		}
	}

	// Calculate alphas of the mesh
	MatrixXf alphasMeshes = eigen_vectors.transpose() * centS;

	// Calculate matrix M
	JacobiSVD<MatrixXf> svd(features, ComputeFullU | ComputeFullV);
	MatrixXf U = svd.matrixU();
	MatrixXf V = svd.matrixV();
	VectorXf singularValues = svd.singularValues();
	MatrixXf M_SigmaInv(MatrixXf::Zero(_n_meshes, n_controllers + 1));
	for (int i = 0; i < min(n_controllers + 1, _n_meshes); i++) {
		if (abs(singularValues(i)) < 1e-06) {
			M_SigmaInv(i, i) = 0;
		}
		else {
			M_SigmaInv(i, i) = 1 / singularValues(i);
		}
	}

	MatrixXf featuresInv = V * M_SigmaInv * U.transpose();
	MatrixXf M_feature2Alpha = alphasMeshes * featuresInv;

	// Calculate initial features
	VectorXf initialFeatures = features.col(0);

	ofstream out(_feature_filename_url, ios::out | ios::trunc);

	out.write((char*)(&n_controllers), sizeof(int));
	out.write((char*)M_feature2Alpha.data(),
		degrees_freedom*(n_controllers + 1)*sizeof(MatrixXf::Scalar));
	out.write((char*)initialFeatures.data(),
		(n_controllers + 1)*sizeof(MatrixXf::Scalar));
	out.close();

}
//=============================================================================
void PCA::centerModel()
{
	float max_x = numeric_limits<float>::min(),
		max_y = numeric_limits<float>::min(),
		max_z = numeric_limits<float>::min();
	float min_x = numeric_limits<float>::max(),
		min_y = numeric_limits<float>::max(),
		min_z = numeric_limits<float>::max();
	for (int i = 0; i < mean_model.size(); i += 3)
	{
		max_x = max(max_x, mean_model(i));
		max_y = max(max_y, mean_model(i + 1));
		max_z = max(max_z, mean_model(i + 2));
		min_x = min(min_x, mean_model(i));
		min_y = min(min_y, mean_model(i + 1));
		min_z = min(min_z, mean_model(i + 2));
	}

	Vector3f center;
	center(0) = min_x + (max_x - min_x) / 2;
	center(1) = min_y + (max_y - min_y) / 2;
	center(2) = min_z + (max_z - min_z) / 2;

	for (int i = 0; i < mean_model.size(); i += 3)
	{
		mean_model(i) = mean_model(i) - center(i % 3);
		mean_model(i + 1) = mean_model(i + 1) - center((i + 1) % 3);
		mean_model(i + 2) = mean_model(i + 2) - center((i + 2) % 3);
		first_model(i) = first_model(i) - center(i % 3);
		first_model(i + 1) = first_model(i + 1) - center((i + 1) % 3);
		first_model(i + 2) = first_model(i + 2) - center((i + 2) % 3);
	}
}
//=============================================================================