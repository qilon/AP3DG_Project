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
#define UPPER_F 70
#define LOWER_F 102
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
using namespace std;
//=============================================================================
const static char* WINDOW_TITLE = "AP3DG";
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 800;
//=============================================================================
class GLViewer
{
public:
	static int moving, beginx, beginy;
	static GLfloat eye[3];
	static GLfloat center[3];
	static GLfloat fieldView;
	static GLfloat aspectRatio;
	static GLfloat depthNear, depthFar;

	static GLfloat anglex;   /* in degrees */
	static GLfloat angley;   /* in degrees */

	static GLfloat light_diffuse[4];  /* Red diffuse light. */
	static GLfloat light_ambient[4];  /* Grey ambient light. */
	static GLfloat light_position[4];  /* Infinite light location. */

	static GLfloat background_colour[4];  /* Background colour. */

	static PolyMesh mesh;

	//=========================================================================
	static void initialize(int *argcp, char **argv);

	static void setMesh(PolyMesh& _mesh);

	static void drawModel(void);
	static void display(void);
	static void mouse(int button, int state, int x, int y);
	static void motion(int x, int y);
	static void init(void);

	static void Key(unsigned char key, int x, int y);

	static void run();
};
//=============================================================================
