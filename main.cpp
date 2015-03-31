/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees
and is provided without guarantee or warrantee expressed or
implied. This program is -not- in the public domain. */

/* Modified by Jan Kautz to match the Coursework */

/* MACOS compile: g++ cube.c -framework OpenGL -framework glut */

#include "GLViewer.h"

using namespace std;
//=============================================================================
const static bool READ_VERTEX_COLORS = false;

const static char* PLY_FILENAME = "./_models/mesh000.ply";
//=============================================================================
int main(int argc, char **argv)
{
	PolyMesh mesh;
	loadMesh(mesh, PLY_FILENAME, READ_VERTEX_COLORS);

	GLViewer viewer;
	viewer.initialize(&argc, argv);
	viewer.setMesh(mesh);
	viewer.run();

	return 0;
}
//=============================================================================