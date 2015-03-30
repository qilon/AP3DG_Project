#include <stdlib.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GL/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "PolyMesh.h"
//=============================================================================
#define UPPER_A 65
#define LOWER_A 97
#define UPPER_D 68
#define LOWER_D 100
#define UPPER_E 69
#define LOWER_E 101
#define UPPER_M 77
#define LOWER_M 109
#define UPPER_Q 81
#define LOWER_Q 113
#define UPPER_R 82
#define LOWER_R 114
#define UPPER_S 83
#define LOWER_S 115
#define UPPER_T 84
#define LOWER_T 116
#define UPPER_W 87
#define LOWER_W 119
//=============================================================================
const static bool READ_VERTEX_COLORS = false;
const static bool WRITE_RESULT = true;

static char* PLY_FILENAME = "./_models/sphere.ply";

using namespace std;
//=============================================================================
int moving, beginx, beginy;
GLfloat anglex = 0;   /* in degrees */
GLfloat angley = 0;   /* in degrees */

GLfloat light_diffuse[] = { 1.0, 0.0, 0.0, 1.0 };  /* Red diffuse light. */
GLfloat light_ambient[] = { 0.1, 0.1, 0.1, 1.0 };  /* Grey ambient light. */
GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };  /* Infinite light location. */

GLfloat background_colour[] = { 0.5f, 0.5f, 0.5f, 0.0f };  /* Background colour. */

GLfloat eye[3] = { 0.0f, 0.0f, 5.0f };
GLfloat center[3] = { 0.0f, 0.0f, 0.0f };

PolyMesh mesh;
//=============================================================================
void initViewer(int *argcp, char **argv);
//=============================================================================