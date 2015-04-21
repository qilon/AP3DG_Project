#include "MyMesh.h"
//=============================================================================
int readMesh(MyMesh& mesh, const char* filename, bool read_vertex_colors)
{
	OpenMesh::IO::Options ropt, wopt;

	if (read_vertex_colors)
	{
		mesh.request_vertex_colors();
		ropt += OpenMesh::IO::Options::VertexColor;
	}

	std::cout << "Reading file... " << std::endl;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "Could not read file: " << filename << std::endl << std::endl;

		return -1;
	}

	return 0;
}
//=============================================================================
Eigen::MatrixXf mesh2EigenMatrix(const MyMesh& mesh) {
	Eigen::MatrixXf cloud(3, mesh.n_vertices());
	// Iterator at the beginning of the vertices list:
	MyMesh::ConstVertexIter It1(mesh.vertices_begin());
	// Iterator at the end of the vertices list:
	MyMesh::ConstVertexIter It2(mesh.vertices_end());
	// For loop that goes through all the vertices in the mesh:
	for (size_t i = 0; It1 != It2; It1++, i++)
	{
		Vec3f vx = vector_cast<Vec3f>(mesh.point(*It1));
		// Each vector of the mesh is added to the cloud:
		cloud(0, i) = vx[0];
		cloud(1, i) = vx[1];
		cloud(2, i) = vx[2];
	}
	return cloud;
}
//=============================================================================
MyMesh eigenMatrix2Mesh(const MyMesh& original, const MatrixXf& inputMatrix) {

	MyMesh mesh = original;
	int nPts = inputMatrix.cols();

	MyMesh::ConstVertexIter vIt1(mesh.vertices_begin());
	MyMesh::ConstVertexIter vIt2(mesh.vertices_end());

	for (size_t i = 0; vIt1 != vIt2; vIt1++, i++) {
		Eigen::Vector3f point = inputMatrix.col(i);
		mesh.point(*vIt1)[0] = point[0];
		mesh.point(*vIt1)[1] = point[1];
		mesh.point(*vIt1)[2] = point[2];
	}
	return mesh;
}
//==============================================================================
void addNoise(string _models_url_preffix, string _models_url_suffix, int first_index, int _n_meshes, string _models_noise_url)
{
	MyMesh mesh;
	for (int iMesh = 0; iMesh < _n_meshes; iMesh++){
		// Load mesh
		string index;
		stringstream convert;
		convert << first_index + iMesh;
		index = convert.str();
		cout << "Reading mesh #" + index << endl;
		if (!OpenMesh::IO::read_mesh(mesh, _models_url_preffix + index + _models_url_suffix))
		{
			std::cerr << "Cannot read mesh #" + index << std::endl;
		}

		MatrixXf model = mesh2EigenMatrix(mesh);
		float mean = 0.0;
		float sigma = 0.005;
		default_random_engine generator;
		normal_distribution<float> gauss(mean, sigma);
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = 0; j < model.cols(); j++) {
				model(i, j) += gauss(generator);
			}
		}
		mesh = eigenMatrix2Mesh(mesh, model);

		cout << "Writing mesh #" + index << endl;
		if (!OpenMesh::IO::write_mesh(mesh, _models_noise_url + index + _models_url_suffix))
		{
			std::cerr << "Cannot write mesh #" + index << std::endl;
		}
	}
}