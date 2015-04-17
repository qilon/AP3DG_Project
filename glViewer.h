#include <stdlib.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GL/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/glui.h>
#endif

#include "pca.h"
#include "Feature.h"
//=============================================================================
#define UPPER_A 65
#define LOWER_A 97
#define UPPER_C 67
#define LOWER_C 99
#define UPPER_D 68
#define LOWER_D 100
#define UPPER_E 69
#define LOWER_E 101
#define UPPER_F 70
#define LOWER_F 102
#define UPPER_I 73
#define LOWER_I 105
#define UPPER_K 75
#define LOWER_K 107
#define UPPER_L 76
#define LOWER_L 108
#define UPPER_M 77
#define LOWER_M 109
#define UPPER_O 79
#define LOWER_O 111
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
using namespace Eigen;
//=============================================================================
class GLViewer
{
private:
	const static char* WINDOW_TITLE;
	const static int WINDOW_WIDTH;
	const static int WINDOW_HEIGHT;

	const static GLfloat FRUSTUM_COEFF;
	const static GLfloat DEPTH_NEAR;
	const static GLfloat DEPTH_FAR;

	const static GLfloat ZOOM_INCR;
	
	const static GLfloat EYE_DISTANCE;

	const static GLfloat RADIUS_OFFSET;
	const static GLint CIRCLE_NUM_LINES;
	const static Vector3f CIRCLE_XY_COLOR;
	const static Vector3f CIRCLE_XZ_COLOR;
	const static Vector3f CIRCLE_YZ_COLOR;

	const static GLfloat LIGHT_AMBIENT[4];  /* Grey ambient light. */
	const static GLfloat LIGHT_POSITION[4];  /* Infinite light location. */
	const static GLfloat BACKGROUND_COLOUR[4];  /* Background colour. */
	const static Vector3f GLViewer::MODEL_COLOR;

	static int moving, beginx, beginy;
	static Vector3f eye;
	static Vector3f center;
	static Vector3f up;
	static GLfloat aspectRatio;
	static GLfloat depthNear, depthFar;

	static GLfloat anglex;   /* in degrees */
	static GLfloat angley;   /* in degrees */

	static MyMesh mesh;

	static GLfloat radius;
	static int showCircles;

	static int idxFeature;

	static PCA pca;

	static int nFeatures;
	static float* features;

	static float translation[];
	static float rotation[];

	static int window_id;
	static GLUI* glui;
	static GLUI_Translation* glui_trans;
	static GLUI_Translation* glui_zoom;
	static GLUI_Checkbox* glui_check_circles;

	//=========================================================================

	static void initGLUT(int *argc, char **argv);
	static void initGLUI(void);

	static void initGLUIComponents(void);
	static void initGLUIFeatures(Feature* _features, int _nFeatures);

	static void display(void);
	static void reshape(int x, int y);
	static void mouse(int button, int state, int x, int y);
	static void motion(int x, int y);
	static void key(unsigned char key, int x, int y);
	static void zoom(GLfloat distance);

	static void drawCircle(GLfloat radius, Vector3f center, GLint plane,
		GLint numLines, Vector3f color);
	static void calculateRadius();
	static void idle(void);

	static void drawModel(void);
	static void drawText(const char *text, int length, int x, int y);

public:
	static void initialize(int *argcp, char **argv);
	static void setMesh(MyMesh& _mesh);
	static void run();
	static void loadPCA(string _pca_filename_url, string _features_filename_url);
};
//=============================================================================
