#include <iostream>
using namespace std;

// -------------------- Eigen
// Doc: http://eigen.tuxfamily.org/dox/
#include <Eigen/Dense>

// -------------------- OpenMesh
// Doc: http://www.openmesh.org/media/Documentations/OpenMesh-Doc-Latest/index.html
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
using namespace OpenMesh;

typedef PolyMesh_ArrayKernelT<> MyMesh;

// ----------------------------------------------------------------------------

Eigen::MatrixXf mesh2EigenMatrix(MyMesh mesh);
MyMesh eigenMatrix2Mesh(MyMesh original, Eigen::MatrixXf inputMatrix);