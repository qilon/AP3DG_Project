#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
//=============================================================================
using namespace OpenMesh;
using namespace Eigen;
//=============================================================================
struct Traits
{
	/// The default coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3f  Point;

	/// The default normal type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3f  Normal;

	/// The default 1D texture coordinate type is float.
	typedef float  TexCoord1D;
	/// The default 2D texture coordinate type is OpenMesh::Vec2f.
	typedef OpenMesh::Vec2f  TexCoord2D;
	/// The default 3D texture coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3f  TexCoord3D;

	/// The default texture index type
	typedef int TextureIndex;

	/// The default color type is OpenMesh::Vec3uc.
	typedef OpenMesh::Vec4f Color;

#ifndef DOXY_IGNORE_THIS
	VertexTraits{};
	HalfedgeTraits{};
	EdgeTraits{};
	FaceTraits{};
#endif

	VertexAttributes(0);
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
	EdgeAttributes(0);
	FaceAttributes(0);
};
//=============================================================================
typedef OpenMesh::PolyMesh_ArrayKernelT<Traits> MyMesh;
//=============================================================================
int readMesh(MyMesh& mesh, const char* filename, bool read_vertex_colors = false);
//=============================================================================
Eigen::MatrixXf mesh2EigenMatrix(const MyMesh& mesh);
//=============================================================================
MyMesh eigenMatrix2Mesh(const MyMesh& original, const MatrixXf& inputMatrix);
//=============================================================================
