#include <stdlib.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GL/glew.h>
#include <GLUT/glut.h>
#include <GL/glui.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/glui.h>
#endif

#include "pca.h"
#include <queue>
//=============================================================================
/**** KEYBOARD KEYS ****/
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

/* MODES */
#define GENERATE_MODE		0
#define RECONSTRUCT_MODE	1

/* TYPES OF POINTS */
#define UNKNOWN_POINT	0
#define KNOWN_POINT		1
//=============================================================================
using namespace std;
using namespace Eigen;
//=============================================================================
class GLViewer
{
private:
	/**** CONSTANTS ****/

	/* WINDOW PARAMETERS */
	const static char* WINDOW_TITLE;
	const static int WINDOW_WIDTH;
	const static int WINDOW_HEIGHT;

	/* PERSPECTIVE PROPERTIES */
	const static GLfloat FRUSTUM_COEFF;
	const static GLfloat DEPTH_NEAR;
	const static GLfloat DEPTH_FAR;
	const static GLfloat EYE_DISTANCE; /* initial eye z-position  */
	const static GLfloat ZOOM_INCR;	/* zoom increment on */

	/* VIEWER COLOURS */
	const static GLfloat LIGHT_AMBIENT[4];
	const static GLfloat LIGHT_POSITION[4];
	const static GLfloat BACKGROUND_COLOUR[4];
	const static MyMesh::Color MODEL_COLOR;

	/* GUIDANCE CIRCLES PARAMETERS */
	const static GLfloat RADIUS_OFFSET;
	const static GLint CIRCLE_NUM_LINES;
	const static GLfloat CIRCLE_XY_COLOR[3];
	const static GLfloat CIRCLE_XZ_COLOR[3];
	const static GLfloat CIRCLE_YZ_COLOR[3];

	/* GLUI CONTROL PARAMETERS */
	const static float TRANSLATION_SPEED;
	const static float ZOOM_SPEED;
	const static float ROTATION_SPIN_FACTOR;

	/* MODE BUTTON TEXT */
	const static char* GENERATE_MODE_TEXT;
	const static char* RECONSTRUCT_MODE_TEXT;

	/* MESH RECONSTRUCTION */
	const static MyMesh::Color RECONSTRUCTED_POINT_COLOR;
	const static MyMesh::Color SELECTED_INDEX_COLOR;
	const static int REMOVE_VERTEX_INDEX;
	const static int REMOVE_N_RINGS;
	const static int REMOVE_MAX_RINGS;

	//=========================================================================

	/**** VARIABLES ****/

	/* MAIN VARIABLES */
	static MyMesh mesh; /* 3d mesh */
	static PCA pca; /* PCA model */
	static float* features; /* feature values */
	static int mode;
	static MyMesh recons_mesh;
	static VectorXi points_state;
	static int remove_vertex_index;
	static int remove_n_rings;

	/* VIEW VARIABLES */
	static GLfloat eye[3]; /* eye position*/
	static GLfloat aspectRatio; /* view aspect ratio*/

	/* ROTATION AND TRANSLATION MATRIXES*/
	static float translation[];
	static float rotation[];

	/* MOUSE CONTROL VARIABLES */
	static int moving;
	static int beginx, beginy;
	static float translationRatio;

	/* GUIDANCE CIRCLES VARIABLES */
	static GLfloat radius; /* radius of circles */
	static int showCircles; /* show or hide circles */

	/* GLUI COMPONENTS */
	static int window_id;
	static GLUI* glui;
	static GLUI_Translation* glui_trans;
	static GLUI_Translation* glui_zoom;
	static GLUI_Checkbox* glui_check_circles;
	static GLUI_Rollout* features_rollout;
	static GLUI_Rollout* recons_rollout;
	static GLUI_Button* glui_modeButton;

	//=========================================================================

	/**** PRIVATE FUNCTIONS ****/

	/* INITIALIZATION METHODS */
	static void initGLUT(int *argc, char **argv);
	static void initGLUI(void);
	static void initGLUIComponents(void);
	static void initGLUIFeatures(FeatureConfig* _features, int _nFeatures);
	static void initGLUIReconstruction();
	static void initGLUIControlPanel();

	/* GLUT AND GLUI FUNCTIONS */
	static void display(void);
	static void reshape(int x, int y);
	static void mouse(int button, int state, int x, int y);
	static void motion(int x, int y);
	static void key(unsigned char key, int x, int y);
	static void idle(void);
	static void modeButtonCallback(int state);
	static void removePointsButtonCallback(int state);
	static void reconstructButtonCallback(int state);

	/* DRAWING FUNCTIONS */
	static void drawCircle(GLfloat _radius, GLint _plane, GLint _numLines,
		const GLfloat* _color);
	static void drawModel(void);

	/* OTHER FUNCTIONS */
	static void zoom(GLfloat distance); /* increase or decrease eye depth */
	static void calculateRadius(); /* updates circles radius based on mesh size */
	static void updateFeature(int _idxFeature); /* update feature in pca and in mesh */
	static void updateMode();
	static void removeReconsMeshRegion(int _vertex_idx, int _n_rings);
	static void reconstruct();

public:
	static void initialize(int *argcp, char **argv);
	static void loadMesh(string _mesh_filename);
	static void loadPCA(string _pca_filename, string _features_filename);
	static void run();
};
//=============================================================================
