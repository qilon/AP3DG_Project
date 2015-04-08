#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "convert.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace std;

void pca(int nMeshes, string ply_models_url_preffix, string pca_result_url);