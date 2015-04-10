/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees
and is provided without guarantee or warrantee expressed or
implied. This program is -not- in the public domain. */

/* Modified by Jan Kautz to match the Coursework */

/* MACOS compile: g++ cube.c -framework OpenGL -framework glut */

#include "GLViewer.h"
//=============================================================================
// PCA PARAMETERS
const static int NUM_MESHES = 71;
const static string PLY_MODELS_URL_PREFFIX = "_models/scapecomp/mesh";
const static string PCA_RESULT_URL = "./_data/pca_result.dat";

// VIEWER TEST PARAMETERS
const static bool READ_VERTEX_COLORS = false;
const static char* PLY_FILENAME = "./_models/scapecomp/mesh0.ply";

//=============================================================================
int main(int argc, char **argv)
{
	
	/* Generate the eigen vectors and eigen values of the covariance matrix for 
	the 71 meshes and keeps the 70 highest ones */
	//PCA pca = PCA(NUM_MESHES, PLY_MODELS_URL_PREFFIX);
	//pca.write(PCA_RESULT_URL);
	

	/* Mesh viewer*/
	GLViewer viewer;
	viewer.initialize(&argc, argv);

	// Reads a mesh to get the vertices connections
	MyMesh mesh;
	loadMesh(mesh, PLY_FILENAME, READ_VERTEX_COLORS);
	viewer.setMesh(mesh);

	// Loads PCA info
	viewer.loadPCA(PCA_RESULT_URL);

	viewer.run();
	

	return 0;
}
//=============================================================================