/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees
and is provided without guarantee or warrantee expressed or
implied. This program is -not- in the public domain. */

/* Modified by Jan Kautz to match the Coursework */

/* MACOS compile: g++ cube.c -framework OpenGL -framework glut */

#include "GLViewer.h"
#include "pca.h"

using namespace std;
//=============================================================================
// PCA PARAMETERS
const static int NUM_MESHES = 3;
const static string PLY_MODELS_URL_PREFFIX = "_models/scapecomp/mesh";
const static string PCA_RESULT_URL = "pca_result.txt";

// VIEWER TEST PARAMETERS
const static bool READ_VERTEX_COLORS = false;
const static char* PLY_FILENAME = "./_models/mesh000.ply";
//=============================================================================
int main(int argc, char **argv)
{
	// pca(NUM_MESHES, PLY_MODELS_URL_PREFFIX, PCA_RESULT_URL);
	
	PolyMesh mesh;
	loadMesh(mesh, PLY_FILENAME, READ_VERTEX_COLORS);

	GLViewer viewer;
	viewer.initialize(&argc, argv);
	viewer.setMesh(mesh);
	viewer.run();

	return 0;
}
//=============================================================================