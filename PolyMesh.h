#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <random>
//=============================================================================
#define UL 1
#define LB 2
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
	typedef OpenMesh::Vec3uc Color;

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
typedef OpenMesh::PolyMesh_ArrayKernelT<Traits> PolyMesh;
//=============================================================================
Eigen::ArrayXf meanCurvature(PolyMesh* mesh, PolyMesh* colored_mesh,
	int laplacian_type);
//=============================================================================
Eigen::ArrayXf gaussianCurvature(PolyMesh* mesh, PolyMesh* colored_mesh);
//=============================================================================
void principalK1K2(PolyMesh* mesh, Eigen::ArrayXf H,
	Eigen::ArrayXf K, Eigen::ArrayXf *k1, Eigen::ArrayXf *k2,
	PolyMesh* k1_mesh, PolyMesh* k2_mesh);
//=============================================================================
PolyMesh explicitSmoothing(PolyMesh* mesh, float diffusion_const,
	int num_iter, int laplacian_type, PolyMesh::Color color);
//=============================================================================
PolyMesh implicitSmoothing(PolyMesh* mesh, float diffusion_const,
	int num_iter, PolyMesh::Color color);
//=============================================================================
void addNoise(PolyMesh* mesh, float mean, float sigma, PolyMesh::Color color);
//=============================================================================