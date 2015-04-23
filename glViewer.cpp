#include "GLViewer.h"
//=============================================================================
/**** CONSTANTS****/

/* WINDOW PARAMETERS */
const char* GLViewer::WINDOW_TITLE = "AP3DG - Qi Liu & Pierre Jestin";
const int GLViewer::WINDOW_WIDTH = 1000;
const int GLViewer::WINDOW_HEIGHT = 800;

/* PERSPECTIVE PROPERTIES */
const GLfloat GLViewer::FRUSTUM_COEFF = 0.025f;
const GLfloat GLViewer::DEPTH_NEAR = 0.1f;
const GLfloat GLViewer::DEPTH_FAR = 50.0f;
const GLfloat GLViewer::EYE_DISTANCE = 5.0f;
const GLfloat GLViewer::ZOOM_INCR = 0.2f;

/* VIEWER COLOURS */
const GLfloat GLViewer::LIGHT_AMBIENT[4] = { 0.1f, 0.1f, 0.1f, 1.f };
const GLfloat GLViewer::LIGHT_POSITION[4] = { .5f, .5f, 1.f, 0.f };
const GLfloat GLViewer::BACKGROUND_COLOUR[4] = { 0.9f, 0.9f, 0.9f, 1.f };
const GLfloat GLViewer::MODEL_COLOR[3] = { .9f, .79f, .69f };

/* GUIDANCE CIRCLES PARAMETERS */
const GLfloat GLViewer::RADIUS_OFFSET = .1f;
const GLint GLViewer::CIRCLE_NUM_LINES = 100;
const GLfloat GLViewer::CIRCLE_XY_COLOR[3] = { 1.f, 0.f, 0.f };
const GLfloat GLViewer::CIRCLE_XZ_COLOR[3] = { 0.f, 1.f, 0.f };
const GLfloat GLViewer::CIRCLE_YZ_COLOR[3] = { 0.f, 0.f, 1.f };

/* GLUI CONTROL PARAMETERS */
const float GLViewer::TRANSLATION_SPEED = .005f;
const float GLViewer::ZOOM_SPEED = .005f;
const float GLViewer::ROTATION_SPIN_FACTOR = .98f;
//=============================================================================
/**** VARIABLES ****/

/* MAIN VARIABLES */
PCA GLViewer::pca = PCA();
MyMesh GLViewer::mesh;
float* GLViewer::features;
float* GLViewer::alphas;

/* VIEW VARIABLES */
GLfloat GLViewer::eye[3] = { 0.f, 0.f, EYE_DISTANCE };
GLfloat GLViewer::aspectRatio = 1.f;

/* ROTATION AND TRANSLATION MATRIXES*/
float GLViewer::translation[3] = { 0.f, 0.f, 0.f };
float GLViewer::rotation[16] = {
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 1
};

/* MOUSE CONTROL VARIABLES */
int GLViewer::moving = 0;
int GLViewer::beginx = 0;
int GLViewer::beginy = 0;
float GLViewer::translationRatio = 0.0025f; // TODO: translation should depend on eye z position

/* GUIDANCE CIRCLES VARIABLES */
GLfloat GLViewer::radius = 1.f;
int GLViewer::showCircles = 1;

/* GLUI COMPONENTS */
GLUI* GLViewer::glui;
int GLViewer::window_id;
GLUI_Translation* GLViewer::glui_trans;
GLUI_Translation* GLViewer::glui_zoom;
GLUI_Checkbox* GLViewer::glui_check_circles;
//=============================================================================
void GLViewer::initGLUT(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - WINDOW_WIDTH) / 2,
		10);

	window_id = glutCreateWindow(WINDOW_TITLE);

	glutDisplayFunc(display);
	glutMotionFunc(motion);

	glewInit();

	/* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_POSITION, LIGHT_POSITION);
	glEnable(GL_COLOR_MATERIAL);
	glLightfv(GL_LIGHT0, GL_AMBIENT, LIGHT_AMBIENT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glClearColor(BACKGROUND_COLOUR[0], BACKGROUND_COLOUR[1], 
		BACKGROUND_COLOUR[2], BACKGROUND_COLOUR[3]);

	/* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);

	/* Anti-aliasing*/
	glEnable(GL_MULTISAMPLE);
}
//=============================================================================
void GLViewer::initGLUI(void)
{
	glui = GLUI_Master.create_glui_subwindow(window_id,
		GLUI_SUBWINDOW_RIGHT);
	glui->set_main_gfx_window(window_id);

	GLUI_Master.set_glutReshapeFunc(reshape);
	GLUI_Master.set_glutKeyboardFunc(key);
	GLUI_Master.set_glutSpecialFunc(NULL);
	GLUI_Master.set_glutMouseFunc(mouse);
	GLUI_Master.set_glutIdleFunc(idle);
}
//=============================================================================
/*
 * Function that initialize the GLUI components
 * Adds the following elements:
 * - list of feature spinners from the info in pca
 * - translation tool
 * - zoom tool
 * - rotation tool
 * - check for showing or hiding guidance circles
 */
void GLViewer::initGLUIComponents(void)
{
	/* Load features from PCA */
	int n_features = pca.getNFeatures();
	initGLUIFeatures(pca.getInitialFeatures(), n_features);
	initGLUIAlphas(pca.getAlphas());

	/* Control Panel */
	GLUI_Panel* control_panel = glui->add_panel("Controls", GLUI_PANEL_NONE);

	GLUI_Panel* trans_panel = glui->add_panel_to_panel(control_panel, "Translations",
		GLUI_PANEL_NONE);

	/* Translation */
	glui_trans =
		new GLUI_Translation(trans_panel, "Translate", GLUI_TRANSLATION_XY, translation);
	glui_trans->set_speed(TRANSLATION_SPEED);

	glui->add_column_to_panel(trans_panel, 0);

	/* Zoom */
	glui_zoom =
		new GLUI_Translation(trans_panel, "Zoom", GLUI_TRANSLATION_Z, &translation[2]);
	glui_zoom->set_speed(ZOOM_SPEED);

	/* Rotation */
	GLUI_Rotation *glui_rot = new GLUI_Rotation(control_panel, "Rotate", rotation);
	glui_rot->set_spin(ROTATION_SPIN_FACTOR);

	/* Circles check */
	glui_check_circles =
		new GLUI_Checkbox(control_panel, "Guidance circles", &showCircles);
	glui_check_circles->set_alignment(GLUI_ALIGN_RIGHT);
}
//=============================================================================
/*
 * Function that creates a spinner for each of the features defined in pca
 */
void GLViewer::initGLUIFeatures(FeatureConfig* _features, int _nFeatures)
{
	GLUI_Panel *features_panel = new GLUI_Rollout(glui, "Features", true);
	features = new float[_nFeatures];
	for (int i = 0; i < _nFeatures; i++)
	{
		features[i] = _features[i].init_value;
		GLUI_Spinner *spinner =
			new GLUI_Spinner(features_panel, _features[i].name,
			&features[i], i, updateFeature);
		spinner->set_float_limits(_features[i].min_value,
			_features[i].max_value);
		spinner->set_float_val(_features[i].init_value);
		//spinner->set_speed(_features[i].getIncrValue());
		spinner->set_alignment(GLUI_ALIGN_RIGHT);
	}
}
//=============================================================================
void GLViewer::initGLUIAlphas(VectorXf _alphas)
{
	GLUI_Panel *alphas_panel = new GLUI_Rollout(glui, "Alphas", false);
	int degrees_freedom = _alphas.size();
	alphas = new float[degrees_freedom];
	for (int i = 0; i < NUM_ALPHAS_DISPLAYED; i++)
	{
		alphas[i] = _alphas(i);
		stringstream ss;
		ss << "Alpha No." << i;
		GLUI_Spinner *spinner =
			new GLUI_Spinner(alphas_panel, "Alpha No.",//ss.str(),
			&alphas[i], i, updateAlpha);
		spinner->set_float_limits(ALPHA_MIN,
			ALPHA_MAX);
		spinner->set_float_val(_alphas(i));
		//spinner->set_speed(_features[i].getIncrValue());
		spinner->set_alignment(GLUI_ALIGN_RIGHT);
	}
}
//=============================================================================
void GLViewer::display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-aspectRatio*FRUSTUM_COEFF, aspectRatio*FRUSTUM_COEFF, 
		-FRUSTUM_COEFF, FRUSTUM_COEFF, DEPTH_NEAR, DEPTH_FAR);

	glMatrixMode(GL_MODELVIEW);

	/* Adjust cube position according to user input */
	/* Note, this is very basic, and suffers from Gimball lock */
	glLoadIdentity();

	gluLookAt(eye[0], eye[1], eye[2],	/* eye */
		0.f, 0.f, 0.f,	/* center */
		0.f, 1.f, 0.f);	/* up is in positive Y direction */

	// Translation
	glTranslatef(translation[0], translation[1], translation[2]);

	// GLUI Rotation
	glMultMatrixf(rotation);

	drawModel();

	// Guidance Circles
	if (showCircles)
	{
		glDisable(GL_LIGHTING);
		drawCircle(radius, 0, CIRCLE_NUM_LINES, CIRCLE_XY_COLOR);
		drawCircle(radius, 1, CIRCLE_NUM_LINES, CIRCLE_XZ_COLOR);
		drawCircle(radius, 2, CIRCLE_NUM_LINES, CIRCLE_YZ_COLOR);
		glEnable(GL_LIGHTING);
	}

	glutSwapBuffers();
}
//=============================================================================
void GLViewer::reshape(int x, int y)
{
	int tx, ty, tw, th;
	GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
	glViewport(tx, ty, tw, th);

	aspectRatio = (float)tw / (float)th;

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::mouse(int button, int state, int x, int y)
{
	// Left button pressed
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		moving = 1;
		beginx = x;
		beginy = y;
	}
	// Left button released
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
		moving = 0;
	}
	// Wheel scroll down
	if (button == 3 && state == GLUT_DOWN) {
		zoom(ZOOM_INCR);
	}
	// Wheel scroll up
	if (button == 4 && state == GLUT_DOWN) {
		zoom(-ZOOM_INCR);
	}
}
//=============================================================================
void GLViewer::motion(int x, int y)
{
	if (moving) {
		float disp_x = (x - beginx) * translationRatio;
		float disp_y = -(y - beginy) * translationRatio;
		beginx = x;
		beginy = y;

		translation[0] += disp_x;
		translation[1] += disp_y;

		glui_trans->set_float_array_val(translation);

		glutPostRedisplay();
	}
}
//=============================================================================
/*
* Controls:
* R/F -> Zoom in/out
* C   -> Show/hide guidance circles
*/
void GLViewer::key(unsigned char key, int x, int y) {
	switch (key) {
	case UPPER_R:
		zoom(ZOOM_INCR);
		break;
	case LOWER_R:
		zoom(ZOOM_INCR);
		break;
	case UPPER_F:
		zoom(-ZOOM_INCR);
		break;
	case LOWER_F:
		zoom(-ZOOM_INCR);
		break;
	case UPPER_C:
		showCircles = !showCircles;
		glui_check_circles->set_int_val(showCircles);
		break;
	case LOWER_C:
		showCircles = !showCircles;
		glui_check_circles->set_int_val(showCircles);
		break;
	default:
		break;
	}
	glutPostRedisplay();
}
//=============================================================================
void GLViewer::idle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != window_id)
		glutSetWindow(window_id);

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::drawCircle(GLfloat _radius, GLint _plane, GLint _numLines,
	const GLfloat* _color)
{
	glPushMatrix();

	// Rotation
	switch (_plane)
	{
	case 0: // xy-plane
		break;
	case 1: // xz-plane
		glRotatef(90.f, 1.f, 0.f, 0.f);
		break;
	case 2: // yz-plane
		glRotatef(90.f, 0.f, 1.f, 0.f);
		break;
	default:
		break;
	}

	// Lines
	glBegin(GL_LINE_LOOP);

	GLfloat x, y, z;
	z = 0.f;
	for (int i = 0; i < _numLines; i++)
	{
		GLfloat angle = 2 * M_PI * i / _numLines;
		x = radius * cos(angle);
		y = radius * sin(angle);
		// Color
		glColor3f(_color[0], _color[1], _color[2]);
		glVertex3f(x, y, z);
	}

	glEnd();

	glPopMatrix();
}
//=============================================================================
void GLViewer::drawModel(void)
{
	glBegin(GL_TRIANGLES);
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); 
		f_it != mesh.faces_end(); ++f_it)
	{
		MyMesh::FaceVertexIter fv_it;
		for (fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			MyMesh::Point p = mesh.point(*fv_it);
			float point[3] {p[0], p[1], p[2]};

			MyMesh::Normal n = mesh.normal(*fv_it);
			float normal[3] {n[0], n[1], n[2]};

			glColor3f(MODEL_COLOR[0], MODEL_COLOR[1], MODEL_COLOR[2]);
			glNormal3fv(normal);
			glVertex3fv(point);
		}
	}
	glEnd();
}
//=============================================================================
void GLViewer::zoom(GLfloat distance)
{
	translation[2] += distance;

	glui_zoom->set_z(translation[2]);
}
//=============================================================================
void GLViewer::calculateRadius()
{
	GLfloat max_x = numeric_limits<float>::min(),
		max_y = numeric_limits<float>::min(),
		max_z = numeric_limits<float>::min();
	GLfloat min_x = numeric_limits<float>::max(),
		min_y = numeric_limits<float>::max(),
		min_z = numeric_limits<float>::max();
	MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
	int i_v = 0;
	for (v_it = mesh.vertices_begin(); v_it != v_end; v_it++)
	{
		MyMesh::Point p = mesh.point(*v_it);
		max_x = max(max_x, p[0]);
		max_y = max(max_y, p[1]);
		max_z = max(max_z, p[2]);
		min_x = min(min_x, p[0]);
		min_y = min(min_y, p[1]);
		min_z = min(min_z, p[2]);
	}

	radius = max(max_x, max(max_y, max_z)) + RADIUS_OFFSET;
}
//=============================================================================
void GLViewer::updateFeature(int _idxFeature)
{
	pca.editFeature(_idxFeature, features[_idxFeature]);
	MatrixXf alphasVector = pca.getAlphas();
	for (int i = 0; i < alphasVector.size(); i++) {
		alphas[i] = alphasVector(i);
	}
	pca.updateMesh(mesh);
}
//=============================================================================
void GLViewer::updateAlpha(int _idxAlpha)
{
	pca.editAlpha(_idxAlpha, alphas[_idxAlpha]);
	MatrixXf featuresVector = pca.getFeatures();
	for (int i = 0; i < featuresVector.size(); i++) {
		features[i] = featuresVector(i);
	}
	pca.updateMesh(mesh);
}
//=============================================================================
void GLViewer::initialize(int *argc, char **argv)
{
	initGLUT(argc, argv);
	initGLUI();
}
//=============================================================================
/*
 * Reads mesh from file and computes the normals if not provided.
 * This mesh is loaded basicly to have faces
 */
void GLViewer::loadMesh(string _mesh_filename)
{
	readMesh(mesh, _mesh_filename.c_str());

	// Add vertex normals
	mesh.request_vertex_normals();

	// Add face normals
	mesh.request_face_normals();

	// If the file did not provide vertex normals, then calculate them
	if (mesh.has_face_normals() && mesh.has_vertex_normals())
	{
		// let the mesh update the normals
		mesh.update_normals();
	}

	// Update guidance circle radius
	calculateRadius();

	glutPostRedisplay();
}
//=============================================================================
/*
 * Reads pca model and the corresponding features data from files
 */
void GLViewer::loadPCA(string _pca_filename, string _features_filename)
{
	pca.readPCA(_pca_filename);
	pca.readFeatures(_features_filename);

	// Update mesh based on the pca
	pca.updateMesh(mesh);

	// Updates guidance circles radius based on the new mesh
	calculateRadius();

	// Initialize UI components based on the features defined in PCA
	initGLUIComponents();
}
//=============================================================================
void GLViewer::run()
{
	glutMainLoop();
}
//=============================================================================
