#include "pca.h"
//=============================================================================
// Stardard, "empty" constructor
PCA::PCA()
{
	degrees_freedom = 0;
	vector_size = 0;
	n_features = 0;
	initialFeatures = nullptr;
}
//=============================================================================
// Constructor that reads both the PCA info and features
PCA::PCA(string _pca_filename_url, string _features_filename_url) : PCA()
{
	readPCA(_pca_filename_url);
	readFeatures(_features_filename_url);
}
//=============================================================================
// Constructor that reads only the PCA info
PCA::PCA(string _pca_filename_url) : PCA()
{
	readPCA(_pca_filename_url);
}
//=============================================================================
// Constructor that computes and writes PCA info out of the models
PCA::PCA(int _n_meshes, string _ply_models_url_preffix, 
	string _ply_models_url_suffix, int first_index, string _pca_filename_url) : PCA()
{
	// Reading meshes:
	MyMesh *meshes = new MyMesh[_n_meshes];

	for (int iMesh = 0; iMesh < _n_meshes; iMesh++){
		string index;
		stringstream convert;
		convert << first_index + iMesh;
		index = convert.str();
		cout << "Reading mesh #" + index << endl;
		if (!OpenMesh::IO::read_mesh(meshes[iMesh], _ply_models_url_preffix 
			+ index + _ply_models_url_suffix))
		{
			std::cerr << "Cannot read mesh #" + index << std::endl;
		}
	}

	// Computing PCA
	computePCA(meshes, _n_meshes, _pca_filename_url);

	// Initialize coefficients
	initAlphas();

	// Deleting meshes
	delete[] meshes;
	meshes = nullptr;
}
//=============================================================================
PCA::PCA(int _n_meshes, string _ply_models_url_preffix, string _ply_models_url_suffix,
	int first_index, string _pca_filename_url, int _gender, string _features_data_filename_url) : PCA()
{
	MatrixXf featuresMeshes = readFeaturesData(_features_data_filename_url);
	int n_meshes_gender = 0;
	MyMesh *meshes = new MyMesh[_n_meshes];

	// Reading meshes:
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++)
	{
		if (featuresMeshes(0, iMesh) == _gender)
		{
			string index;
			stringstream convert;
			convert << first_index + iMesh;
			index = convert.str();
			cout << "Reading mesh #" + index << endl;
			if (!OpenMesh::IO::read_mesh(meshes[iMesh], _ply_models_url_preffix
				+ index + _ply_models_url_suffix))
			{
				std::cerr << "Cannot read mesh #" + index << std::endl;
			}
			n_meshes_gender++;
		}
	}

	// Computing PCA
	computePCA(meshes, n_meshes_gender, _pca_filename_url);

	// Initialize coefficients
	initAlphas();

	// Deleting meshes
	delete[] meshes;
	meshes = nullptr;
}
//=============================================================================
// Destructor
PCA::~PCA()
{
	if (initialFeatures != nullptr) {
		delete[] initialFeatures;
		initialFeatures = nullptr;
	}
}
//=============================================================================
// Function that computes PCA info. Warning: needs a lot of resources.
void PCA::computePCA(MyMesh* meshes, int _n_meshes, string _pca_filename_url)
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

	writePCA(_pca_filename_url);
}
//=============================================================================
// Function that writes PCA info into a file.
void PCA::writePCA(string _pca_filename_url)
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
// Read PCA info from a file.
void PCA::readPCA(string _pca_filename_url)
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

	// Read the eigenvectors matrix
	inPCA.read((char *)eigen_vectors.data(),
		vector_size*degrees_freedom*sizeof(MatrixXf::Scalar));
	// Read the eigenvalues matrix
	inPCA.read((char *)eigen_values.data(),
		degrees_freedom*sizeof(VectorXf::Scalar));
	// Read the mean model
	inPCA.read((char*)mean_model.data(),
		vector_size*sizeof(VectorXf::Scalar));
	// Read the first model
	inPCA.read((char*)first_model.data(),
		vector_size*sizeof(VectorXf::Scalar));

	inPCA.close();
}
//=============================================================================
// Read the table of feature values from the metadata associated with the Body Data set of meshes
MatrixXf PCA::readFeaturesData(string _features_data_filename_url)
{
	// Load features data (e.g. from the body_data metadata)
	ifstream inFeatures(_features_data_filename_url);

	inFeatures >> n_features;
	cout << "Number of features: " << n_features << endl;

	int n_meshes;
	inFeatures >> n_meshes;
	cout << "Number of meshes: " << n_meshes << endl;

	// Read feature labels
	initialFeatures = new FeatureConfig[n_features];
	string line;
	getline(inFeatures, line);
	for (int iFeature = 0; iFeature < n_features; iFeature++) {
		getline(inFeatures, line);
		strcpy(initialFeatures[iFeature].name, line.c_str());
	}

	// Min values
	for (int iFeature = 0; iFeature < n_features; iFeature++) {
		inFeatures >> initialFeatures[iFeature].min_value;
	}

	// Max values
	for (int iFeature = 0; iFeature < n_features; iFeature++) {
		inFeatures >> initialFeatures[iFeature].max_value;
	}

	// Read feature values
	getline(inFeatures, line);
	MatrixXf featuresMeshes(MatrixXf::Ones(n_features + 1, n_meshes));
	for (int iMesh = 0; iMesh < n_meshes; iMesh++) {
		getline(inFeatures, line);
		stringstream iss(line);
		float x;
		int iFeature = 0;
		while (iss >> x) {
			featuresMeshes(iFeature, iMesh) = x;
			iFeature++;
		}
	}

	// Initial value
	for (int iFeature = 0; iFeature < n_features; iFeature++) {
		initialFeatures[iFeature].init_value = featuresMeshes(iFeature, 0);
	}

	return featuresMeshes;
}
//=============================================================================
// Compute the initial features and conversion matrix for the Body Data set of meshes
void PCA::computeFeatures(string _features_data_filename_url, string _ply_models_url_preffix,
	string _ply_models_url_suffix, int first_index, string _feature_filename_url) {

	// First: read the features table from a file
	MatrixXf featuresMeshes = readFeaturesData(_features_data_filename_url);

	// Load the meshes
	int n_meshes = featuresMeshes.cols();
	MyMesh *meshes = new MyMesh[n_meshes];
	for (int iMesh = 0; iMesh < n_meshes; iMesh++){
		string index;
		stringstream convert;
		convert << first_index + iMesh;
		index = convert.str();
		cout << "Reading mesh #" + index << endl;
		if (!OpenMesh::IO::read_mesh(meshes[iMesh], _ply_models_url_preffix + index + _ply_models_url_suffix))
		{
			std::cerr << "Cannot read mesh #" + index << std::endl;
		}
	}

	// Compute centered vectors with the coordinates of the vertices
	int nVert = mean_model.size() / 3;
	MatrixXf centS(3 * nVert, n_meshes);
	for (int iMesh = 0; iMesh < n_meshes; iMesh++) {
		MatrixXf model = mesh2EigenMatrix(meshes[iMesh]);
		// cout << "min x: " << model.row(0).minCoeff() << endl;
		// cout << "max x: " << model.row(0).maxCoeff() << endl;
		for (int iVert = 0; iVert < nVert; iVert++) {
			centS(3 * iVert, iMesh) = model(0, iVert) - mean_model(3 * iVert);
			centS(3 * iVert + 1, iMesh) = model(1, iVert) - mean_model(3 * iVert + 1);
			centS(3 * iVert + 2, iMesh) = model(2, iVert) - mean_model(3 * iVert + 2);
		}
	}

	// Calculate alphas of the mesh
	MatrixXf alphasMeshes = eigen_vectors.transpose() * centS;

	// Calculate matrix M
	JacobiSVD<MatrixXf> svd(featuresMeshes, ComputeFullU | ComputeFullV);
	MatrixXf U = svd.matrixU();
	MatrixXf V = svd.matrixV();
	VectorXf singularValues = svd.singularValues();
	MatrixXf M_SigmaInv(MatrixXf::Zero(n_meshes, n_features + 1));
	for (int i = 0; i < min(n_features + 1, n_meshes); i++) {
		if (abs(singularValues(i)) > 1e-06) {
			M_SigmaInv(i, i) = 1 / singularValues(i);
		}
	}

	// Calculate conversion matrix
	MatrixXf featuresInv = V * M_SigmaInv * U.transpose();
	M_feature2Alpha = alphasMeshes * featuresInv;

	writeFeatures(_feature_filename_url);
}
//=============================================================================
void PCA::computeFeatures(string _features_data_filename_url, string _ply_models_url_preffix,
	string _ply_models_url_suffix, int first_index, string _feature_filename_url, int _gender) {

	// First: read the features table from a file
	MatrixXf featuresMeshes = readFeaturesData(_features_data_filename_url);
	int n_meshes = featuresMeshes.cols();

	int n_meshes_gender = 0;
	MyMesh *meshes = new MyMesh[n_meshes];

	// Reading meshes:
	for (int iMesh = 0; iMesh < n_meshes; iMesh++)
	{
		if (featuresMeshes(0, iMesh) == _gender)
		{
			string index;
			stringstream convert;
			convert << first_index + iMesh;
			index = convert.str();
			cout << "Reading mesh #" + index << endl;
			if (!OpenMesh::IO::read_mesh(meshes[iMesh], _ply_models_url_preffix
				+ index + _ply_models_url_suffix))
			{
				std::cerr << "Cannot read mesh #" + index << std::endl;
			}
			n_meshes_gender++;
		}
	}

	// Compute centered vectors with the coordinates of the vertices
	int nVert = mean_model.size() / 3;
	MatrixXf centS(3 * nVert, n_meshes_gender);
	for (int iMesh = 0; iMesh < n_meshes_gender; iMesh++) {
		MatrixXf model = mesh2EigenMatrix(meshes[iMesh]);
		for (int iVert = 0; iVert < nVert; iVert++) {
			centS(3 * iVert, iMesh) = model(0, iVert) - mean_model(3 * iVert);
			centS(3 * iVert + 1, iMesh) = model(1, iVert) - mean_model(3 * iVert + 1);
			centS(3 * iVert + 2, iMesh) = model(2, iVert) - mean_model(3 * iVert + 2);
		}
	}

	// Calculate alphas of the mesh
	MatrixXf alphasMeshes = eigen_vectors.transpose() * centS;

	// Calculate matrix M
	JacobiSVD<MatrixXf> svd(featuresMeshes, ComputeFullU | ComputeFullV);
	MatrixXf U = svd.matrixU();
	MatrixXf V = svd.matrixV();
	VectorXf singularValues = svd.singularValues();
	MatrixXf M_SigmaInv(MatrixXf::Zero(n_meshes_gender, n_features + 1));
	for (int i = 0; i < min(n_features + 1, n_meshes_gender); i++) {
		if (abs(singularValues(i)) > 1e-06) {
			M_SigmaInv(i, i) = 1 / singularValues(i);
		}
	}

	// Calculate conversion matrix
	MatrixXf featuresInv = V * M_SigmaInv * U.transpose();
	M_feature2Alpha = alphasMeshes * featuresInv;

	writeFeatures(_feature_filename_url);
}
//=============================================================================
// Compute features corresponding to the Scapecomp set of meshes
void PCA::computeFeatures(int _n_meshes, string _ply_models_url_preffix, string _feature_filename_url) {

	n_features = 8;

	// Load meshes
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

	int nVert = mean_model.size() / 3;
	MatrixXf featuresMeshes(MatrixXf::Ones(n_features + 1, _n_meshes));
	MatrixXf centS(3 * nVert, _n_meshes);
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++){
		MatrixXf model = mesh2EigenMatrix(meshes[iMesh]);

		// Compute the different features from select vertices from the meshes
		// Right arm
		Vector3f uRightArm = model.col(RIGHT_INSIDE_ELBOW) - model.col(RIGHT_AXILLA);
		Vector3f vRightArm = model.col(RIGHT_WAIST) - model.col(RIGHT_AXILLA);
		float angleRightArm = acos(uRightArm.dot(vRightArm) / (uRightArm.norm() * vRightArm.norm()));
		featuresMeshes(0, iMesh) = angleRightArm;

		// Left arm
		Vector3f uLeftArm = model.col(LEFT_INSIDE_ELBOW) - model.col(LEFT_AXILLA);
		Vector3f vLeftArm = model.col(LEFT_WAIST) - model.col(LEFT_AXILLA);
		float angleLeftArm = acos(uLeftArm.dot(vLeftArm) / (uLeftArm.norm() * vLeftArm.norm()));
		featuresMeshes(1, iMesh) = angleLeftArm;

		// Right forearm
		Vector3f uRightElbow = (model.col(RIGHT_SHOULDER) + model.col(RIGHT_AXILLA)
			- model.col(RIGHT_INSIDE_ELBOW) - model.col(RIGHT_ELBOW)) / 2;
		Vector3f vRightElbow = (model.col(RIGHT_INSIDE_ELBOW) + model.col(RIGHT_ELBOW)
			- model.col(RIGHT_BACK_WRIST) - model.col(RIGHT_FRONT_WRIST)) / 2;
		float angleRightElbow = acos(uRightElbow.dot(vRightElbow) / (uRightElbow.norm() * vRightElbow.norm()));
		featuresMeshes(2, iMesh) = angleRightElbow;

		// Left forearm
		Vector3f uLeftElbow = (model.col(LEFT_SHOULDER) + model.col(LEFT_AXILLA)
			- model.col(LEFT_INSIDE_ELBOW) - model.col(LEFT_ELBOW)) / 2;
		Vector3f vLeftElbow = (model.col(LEFT_INSIDE_ELBOW) + model.col(LEFT_ELBOW)
			- model.col(LEFT_BACK_WRIST) - model.col(LEFT_FRONT_WRIST)) / 2;
		float angleLeftElbow = acos(uLeftElbow.dot(vLeftElbow) / (uLeftElbow.norm() * vLeftElbow.norm()));
		featuresMeshes(3, iMesh) = angleLeftElbow;

		// Right leg
		Vector3f uRightLeg = (2 * model.col(RIGHT_WAIST) + model.col(LEFT_WAIST)) / 3
			- (model.col(RIGHT_HIP) + model.col(PERINEUM)) / 2;
		Vector3f vRightLeg = (model.col(RIGHT_HIP) + model.col(PERINEUM)
			- model.col(RIGHT_KNEE) - model.col(RIGHT_INSIDE_KNEE)) / 2;
		float angleRightLeg = acos(uRightLeg.dot(vRightLeg) / (uRightLeg.norm() * vRightLeg.norm()));
		featuresMeshes(4, iMesh) = angleRightLeg;

		// Left leg
		Vector3f uLeftLeg = (model.col(RIGHT_WAIST) + 2 * model.col(LEFT_WAIST)) / 3
			- (model.col(LEFT_HIP) + model.col(PERINEUM)) / 2;
		Vector3f vLeftLeg = (model.col(LEFT_HIP) + model.col(PERINEUM)
			- model.col(LEFT_KNEE) - model.col(LEFT_INSIDE_KNEE)) / 2;
		float angleLeftLeg = acos(uLeftLeg.dot(vLeftLeg) / (uLeftLeg.norm() * vLeftLeg.norm()));
		featuresMeshes(5, iMesh) = angleLeftLeg;

		// Right calf
		Vector3f uRightKnee = (model.col(RIGHT_HIP) + model.col(PERINEUM)
			- model.col(RIGHT_INSIDE_KNEE) - model.col(RIGHT_KNEE)) / 2;
		Vector3f vRightKnee = (model.col(RIGHT_INSIDE_KNEE) + model.col(RIGHT_KNEE)
			- model.col(RIGHT_INSIDE_ANKLE) - model.col(RIGHT_OUTSIDE_ANKLE)) / 2;
		float angleRightKnee = acos(uRightKnee.dot(vRightKnee) / (uRightKnee.norm() * vRightKnee.norm()));
		featuresMeshes(6, iMesh) = angleRightKnee;

		// Left calf
		Vector3f uLeftKnee = (model.col(LEFT_HIP) + model.col(PERINEUM)
			- model.col(LEFT_INSIDE_KNEE) - model.col(LEFT_KNEE)) / 2;
		Vector3f vLeftKnee = (model.col(LEFT_INSIDE_KNEE) + model.col(LEFT_KNEE)
			- model.col(LEFT_INSIDE_ANKLE) - model.col(LEFT_OUTSIDE_ANKLE)) / 2;
		float angleLeftKnee = acos(uLeftKnee.dot(vLeftKnee) / (uLeftKnee.norm() * vLeftKnee.norm()));
		featuresMeshes(7, iMesh) = angleLeftKnee;

		// Compute centered vectors with the coordinates of the vertices
		for (int iVert = 0; iVert < nVert; iVert++) {
			centS(3 * iVert, iMesh) = model(0, iVert) - mean_model(3 * iVert);
			centS(3 * iVert + 1, iMesh) = model(1, iVert) - mean_model(3 * iVert + 1);
			centS(3 * iVert + 2, iMesh) = model(2, iVert) - mean_model(3 * iVert + 2);
		}
	}

	// Calculate initial features
	initialFeatures = new FeatureConfig[n_features];
	initialFeatures[0] = FeatureConfig("Right arm", featuresMeshes(0, 0), -5, 5, 0.1);
	initialFeatures[1] = FeatureConfig("Left arm", featuresMeshes(1, 0), -5, 5, 0.1);
	initialFeatures[2] = FeatureConfig("Right elbow", featuresMeshes(2, 0), -5, 5, 0.1);
	initialFeatures[3] = FeatureConfig("Left elbow", featuresMeshes(3, 0), -5, 5, 0.1);
	initialFeatures[4] = FeatureConfig("Right leg", featuresMeshes(4, 0), -5, 5, 0.1);
	initialFeatures[5] = FeatureConfig("Left leg", featuresMeshes(5, 0), -5, 5, 0.1);
	initialFeatures[6] = FeatureConfig("Right knee", featuresMeshes(6, 0), -5, 5, 0.1);
	initialFeatures[7] = FeatureConfig("Left knee", featuresMeshes(7, 0), -5, 5, 0.1);

	// Calculate alphas of the mesh
	MatrixXf alphasMeshes = eigen_vectors.transpose() * centS;

	// Calculate matrix M
	JacobiSVD<MatrixXf> svd(featuresMeshes, ComputeFullU | ComputeFullV);
	MatrixXf U = svd.matrixU();
	MatrixXf V = svd.matrixV();
	VectorXf singularValues = svd.singularValues();
	MatrixXf M_SigmaInv(MatrixXf::Zero(_n_meshes, n_features + 1));
	for (int i = 0; i < min(n_features + 1, _n_meshes); i++) {
		if (abs(singularValues(i)) > 1e-06) {
			M_SigmaInv(i, i) = 1 / singularValues(i);
		}
	}

	// Calculate conversion matrix
	MatrixXf featuresInv = V * M_SigmaInv * U.transpose();
	M_feature2Alpha = alphasMeshes * featuresInv;

	writeFeatures(_feature_filename_url);
}
//=============================================================================
// Write the conversion matrix and the initial features into a file
void PCA::writeFeatures(string _feature_filename_url)
{
	ofstream out(_feature_filename_url, ios::out | std::ios::binary | ios::trunc);

	// Write the number of features
	out.write((char*)(&n_features), sizeof(int));
	// Write the conversion matrix
	out.write((char*)M_feature2Alpha.data(),
		degrees_freedom*(n_features + 1)*sizeof(MatrixXf::Scalar));
	// Write the initial features under the form of a struct FeatureConfig
	out.write((char*)initialFeatures,
		n_features*sizeof(FeatureConfig));
	out.close();
}
//=============================================================================
// Read the conversion matrix and initial features into a file
void PCA::readFeatures(string _features_filename_url)
{
	// Load features
	ifstream inFeatures(_features_filename_url, ios::in | std::ios::binary);
	
	inFeatures.read((char*)(&n_features), sizeof(int));
	
	cout << "Number of features: " << n_features << endl;
	M_feature2Alpha.resize(degrees_freedom, n_features + 1);
	initialFeatures = new FeatureConfig[n_features];

	// Read conversion matrix
	inFeatures.read((char *)M_feature2Alpha.data(),
		degrees_freedom*(n_features + 1)*sizeof(MatrixXf::Scalar));
	// Read initial features
	inFeatures.read((char *)initialFeatures,
		n_features*sizeof(FeatureConfig));
	inFeatures.close();

	// Initialize the features to their initial values
	initFeatures();
	// Initialize the alphas to those of the first mesh
	initAlphas();
	// Center the mean and first models
	centerModel();
}
//=============================================================================
// Update the positions of the vertices in the mesh from the new features
void PCA::updateMesh(MyMesh& _mesh)
{
	VectorXf new_model = mean_model;

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
	_mesh.update_normals();
}
//=============================================================================
// Center the mean and first model in order for them to be at the center of the window in the viewer
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
		mean_model(i) = mean_model(i) - center(0);
		mean_model(i + 1) = mean_model(i + 1) - center(1);
		mean_model(i + 2) = mean_model(i + 2) - center(2);
		first_model(i) = first_model(i) - center(0);
		first_model(i + 1) = first_model(i + 1) - center(1);
		first_model(i + 2) = first_model(i + 2) - center(2);
	}
}
//=============================================================================
// Give initial values to the alphas, in this case, those of the first mesh
void PCA::initAlphas()
{
	alphas = M_feature2Alpha * features;
}
//=============================================================================
void PCA::editAlpha(int _idxAlpha, float value)
{
	alphas(_idxAlpha) = value;
	features = M_feature2Alpha.fullPivHouseholderQr().solve(alphas);
}
//=============================================================================
VectorXf PCA::getAlphas()
{
	return alphas;
}
//=============================================================================
// Give initial values to the features, from the features file
void PCA::initFeatures()
{
	features.resize(n_features + 1);
	for (int iFeature = 0; iFeature < n_features; iFeature++) {
		features(iFeature) = initialFeatures[iFeature].init_value;
	}
	features(n_features) = 1;
}
//=============================================================================
// Change the value of one feature
void PCA::editFeature(int idxFeature, float new_value)
{
	features(idxFeature) = new_value;
	alphas = M_feature2Alpha * features;
}
//=============================================================================
// Get the value of one feature
VectorXf PCA::getFeatures()
{
	return features;
}
//=============================================================================
// Get the entirety of the initial features
FeatureConfig* PCA::getInitialFeatures()
{
	return initialFeatures;
}
//=============================================================================
// Get the number of features
int PCA::getNFeatures()
{
	return n_features;
}
//=============================================================================
void PCA::updateMesh(MyMesh& _mesh, const VectorXf& _v_mesh, 
	const VectorXi& _points_state, const MyMesh::Color& _color)
{
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());
	int i = 0;
	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it) {
		int v_idx = v_it.handle().idx();
		float x = _v_mesh(3 * v_idx);
		float y = _v_mesh(3 * v_idx + 1);
		float z = _v_mesh(3 * v_idx + 2);
		_mesh.point(*v_it)[0] = x;
		_mesh.point(*v_it)[1] = y;
		_mesh.point(*v_it)[2] = z;
		if (!_points_state(v_idx))
			_mesh.set_color(*v_it, _color);
	}

	/* Update mesh normals */
	_mesh.update_normals();
}
//=============================================================================
void PCA::updateMeshAlphas(MyMesh& _mesh, const VectorXf& _alphas,
	const VectorXi& _points_state, const MyMesh::Color& _color)
{
	VectorXf new_model = mean_model;

	for (int i = 0; i < degrees_freedom; i++)
	{
		new_model += _alphas(i) * eigen_vectors.col(i);
	}

	updateMesh(_mesh, new_model, _points_state, _color);
}
//=============================================================================
/* Function that, given a mesh and a vector identifying which points are given 
and which are unknown, reconstructs the missing points */
void PCA::reconstructMesh(MyMesh& _mesh, const VectorXi& _points_state,
	const MyMesh::Color& _color)
{
	int n_knowns = _points_state.sum(); // Number of known points

	/* Column vector representation of the mesh */
	VectorXf v_mesh = mesh2EigenVector(_mesh);

	/* Vector of the difference of the known points S' of the mesh with the mean 
	model */
	VectorXf diff_s_p(3 * n_knowns);

	// Matrix of known coordinates of eigen vectors - V'
	MatrixXf V_p(3 * n_knowns, degrees_freedom);

	/* Fills the matrix V' and store indexes of unknown points */
	int known_row = 0;
	int unknown_row = 0;
	for (int i = 0; i < _points_state.size(); i++)
	{
		if (_points_state(i))
		{
			V_p.block(3 * known_row, 0, 3, degrees_freedom) = 
				eigen_vectors.block(3 * i, 0, 3, degrees_freedom);
			diff_s_p.segment(3 * known_row, 3) = 
				v_mesh.segment(3 * i, 3) - mean_model.segment(3 * i, 3);
			known_row++;
		}
	}

	/* (V'^T * V' + E^-1)a = V'^T * diff(s') */
	JacobiSVD<MatrixXf> svd(V_p.transpose() * V_p + 
		MatrixXf(VectorXf(eigen_values.array().pow(-2)).asDiagonal()),
		ComputeThinU | ComputeThinV);
	VectorXf alphas = svd.solve(V_p.transpose() * diff_s_p);

	updateMeshAlphas(_mesh, alphas, _points_state, _color);
}
//=============================================================================