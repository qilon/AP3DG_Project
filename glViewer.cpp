#include "GLViewer.h"
//=============================================================================
const char* GLViewer::WINDOW_TITLE = "AP3DG - Qi Liu & Pierre Jestin";
const int GLViewer::WINDOW_WIDTH = 1000;
const int GLViewer::WINDOW_HEIGHT = 800;

const GLfloat GLViewer::FRUSTUM_COEFF = 0.025f;
const GLfloat GLViewer::DEPTH_NEAR = 0.1f;
const GLfloat GLViewer::DEPTH_FAR = 50.0f;

const GLfloat GLViewer::ZOOM_INCR = 0.2f;

const GLfloat GLViewer::EYE_DISTANCE = 5.0f;

const GLfloat GLViewer::RADIUS_OFFSET = .1f;
const GLint GLViewer::CIRCLE_NUM_LINES = 100;
const Vector3f GLViewer::CIRCLE_XY_COLOR = {1.f, 0.f, 0.f};
const Vector3f GLViewer::CIRCLE_XZ_COLOR = { 0.f, 1.f, 0.f };
const Vector3f GLViewer::CIRCLE_YZ_COLOR = { 0.f, 0.f, 1.f };

const GLfloat GLViewer::LIGHT_AMBIENT[4] = { 0.1, 0.1, 0.1, 1.0 };  /* Ambient light. */
const GLfloat GLViewer::LIGHT_POSITION[4] = { .5, .5, 1.0, 0.0 };  /* Infinite light location. */
const GLfloat GLViewer::BACKGROUND_COLOUR[4] = { 0.9f, 0.9f, 0.9f, 1.0f };	/* Background colour. */

const Vector3f GLViewer::MODEL_COLOR = { .9f, .79f, .69f };  /* Diffuse light. */

int GLViewer::moving = 0;
int GLViewer::beginx = 0;
int GLViewer::beginy = 0;

Vector3f GLViewer::eye = { 0.0f, 0.0f, EYE_DISTANCE };
Vector3f GLViewer::center = { 0.0f, 0.0f, 0.0f };
Vector3f GLViewer::up = { 0.0f, 1.0f, 0.0f };

GLfloat GLViewer::aspectRatio = 1.0f;

GLfloat GLViewer::anglex = 0;   /* in degrees */
GLfloat GLViewer::angley = 0;   /* in degrees */

MyMesh GLViewer::mesh;

GLfloat GLViewer::radius = 1.f;
int GLViewer::showCircles = 1;

int GLViewer::idxFeature = 0;

PCA GLViewer::pca = PCA();

GLUI* GLViewer::glui;
int GLViewer::window_id;

float* GLViewer::features;

float GLViewer::translation[3] = { 0.f, 0.f, 0.f };
float GLViewer::rotation[16] = { 
	1, 0, 0, 0, 
	0, 1, 0, 0, 
	0, 0, 1, 0, 
	0, 0, 0, 1 
};

GLUI_Translation* GLViewer::glui_trans;
GLUI_Translation* GLViewer::glui_zoom;
GLUI_Checkbox* GLViewer::glui_check_circles;
//=============================================================================
void GLViewer::initialize(int *argc, char **argv)
{
	initGLUT(argc, argv);
	initGLUI();
}
//=============================================================================
void GLViewer::initGLUT(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - WINDOW_WIDTH) / 2,
		10);
	//glutInitWindowPosition(20, 20);
	window_id = glutCreateWindow(WINDOW_TITLE);

	glutDisplayFunc(display);
	glutMotionFunc(motion);

	glewInit();

	/* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_POSITION, LIGHT_POSITION);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
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

	idxFeature = 0;
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

	gluLookAt(eye[0], eye[1], eye[2],  /* eye */
		center[0], center[1], center[2], /* center */
		up[0], up[1], up[2]);      /* up is in positive Y direction */

	// Translation
	glTranslatef(translation[0], translation[1], translation[2]);

	// Rotation
	glMultMatrixf(rotation);

	drawModel();

	if (showCircles)
	{
		glDisable(GL_LIGHTING);
		drawCircle(radius, center, 0, CIRCLE_NUM_LINES, CIRCLE_XY_COLOR);
		drawCircle(radius, center, 1, CIRCLE_NUM_LINES, CIRCLE_XZ_COLOR);
		drawCircle(radius, center, 2, CIRCLE_NUM_LINES, CIRCLE_YZ_COLOR);
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
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		moving = 1;
		beginx = x;
		beginy = y;
	}
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
		// TODO: change to depend on z translation
		float disp_x = (x - beginx) * 0.0025f;
		float disp_y = -(y - beginy) * 0.0025f;
		beginx = x;
		beginy = y;

		translation[0] += disp_x;
		translation[1] += disp_y;

		glui_trans->set_float_array_val(translation);

		glutPostRedisplay();
	}
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

			glColor3f(MODEL_COLOR(0), MODEL_COLOR(1), MODEL_COLOR(2));
			glNormal3fv(normal);
			glVertex3fv(point);
		}
	}
	glEnd();
}
//=============================================================================
/*
 * Controls:
 * W/S -> Rotate around x axis
 * A/D -> Rotate around y axis
 * Q/E -> Rotate around z axis
 * R/F -> Zoom in/out
 * C   -> Show/hide guidance circles
 * K/L -> Change alpha coefficient
 */
void GLViewer::key(unsigned char key, int x, int y) {
	switch (key) {
	case UPPER_W:
		break;
	case LOWER_W:
		break;
	case UPPER_S:
		break;
	case LOWER_S:
		break;
	case UPPER_A:
		break;
	case LOWER_A:
		break;
	case UPPER_D:
		break;
	case LOWER_D:
		break;
	case UPPER_Q:
		break;
	case LOWER_Q:
		break;
	case UPPER_E:
		break;
	case LOWER_E:
		break;
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
	case UPPER_K:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) - 0.1);
		pca.updateMesh(mesh);
		break;
	case LOWER_K:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) - 0.1);
		pca.updateMesh(mesh);
		break;
	case UPPER_L:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) + 0.1);
		pca.updateMesh(mesh);
		break;
	case LOWER_L:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) + 0.1);
		pca.updateMesh(mesh);
		break;
	case UPPER_I:
		if (idxFeature > 0)
			idxFeature--;
		break;
	case LOWER_I:
		if (idxFeature > 0)
			idxFeature--;
		break;
	case UPPER_O:
		if (idxFeature < pca.getFeatures() - 1) {
			idxFeature++;
		}
		break;
	case LOWER_O:
		if (idxFeature < pca.getFeatures() - 1) {
			idxFeature++;
		}
		break;
	default:
		break;
	}
	glutPostRedisplay();
}
//=============================================================================
void GLViewer::setMesh(MyMesh& _mesh)
{
	mesh = _mesh;

	// Add vertex normals as default property (ref. previous tutorial)
	mesh.request_vertex_normals();

	// Add face normals as default property
	mesh.request_face_normals();

	// If the file did not provide vertex normals, then calculate them
	if (mesh.has_face_normals() && mesh.has_vertex_normals())
	{
		// let the mesh update the normals
		mesh.update_normals();
	}

	calculateRadius();

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::zoom(GLfloat distance)
{
	translation[2] += distance;

	glui_zoom->set_z(translation[2]);
}
//=============================================================================
void GLViewer::run()
{
	glutMainLoop();
}
//=============================================================================
void GLViewer::drawCircle(GLfloat radius, Vector3f center, GLint plane, 
	GLint numLines, Vector3f color)
{
	glPushMatrix();

	// Rotation
	switch (plane)
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
	for (int i = 0; i < numLines; i++)
	{
		GLfloat angle = 2 * M_PI * i / numLines;
		x = radius * cos(angle);
		y = radius * sin(angle);
		// Color
		glColor3f(color(0), color(1), color(2));
		glVertex3f(x, y, z);
	}

	glEnd();

	glPopMatrix();
}
//=============================================================================
void GLViewer::loadPCA(string _pca_filename_url, string _features_filename_url)
{
	pca.readPCA(_pca_filename_url);
	pca.readFeatures(_features_filename_url);

	// Update mesh based on the pca
	pca.updateMesh(mesh);

	calculateRadius();

	// Initialize UI components based on the features defined in PCA
	initGLUIComponents();
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

	radius = max(max_x,max(max_y,max_z)) + RADIUS_OFFSET;
}
//=============================================================================
// https://www.youtube.com/watch?v=elE__Nouv54
void GLViewer::drawText(const char *text, int length, int x, int y) {
	glMatrixMode(GL_PROJECTION);
	double *matrix = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, matrix);
	glLoadIdentity();
	glOrtho(0, 800, 0, 600, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);
	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]);
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(matrix);
	glMatrixMode(GL_MODELVIEW);
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
void GLViewer::initGLUIComponents(void)
{
	/* Load features from PCA */
	int n_features = pca.getFeatures();
	initGLUIFeatures(pca.getInitialFeatures(), n_features);

	/* Control Panel */
	GLUI_Panel* control_panel = glui->add_panel("Controls", GLUI_PANEL_NONE);

	GLUI_Panel* trans_panel = glui->add_panel_to_panel(control_panel, "Translations",
		GLUI_PANEL_NONE);

	glui_trans =
		new GLUI_Translation(trans_panel, "Translate", GLUI_TRANSLATION_XY, translation);
	glui_trans->set_speed(.005);

	glui->add_column_to_panel(trans_panel, 0);

	glui_zoom =
		new GLUI_Translation(trans_panel, "Zoom", GLUI_TRANSLATION_Z, &translation[2]);
	glui_zoom->set_speed(.005);

	GLUI_Rotation *glui_rot = new GLUI_Rotation(control_panel, "Rotate", rotation);
	glui_rot->set_spin(.98);

	glui_check_circles =
		new GLUI_Checkbox(control_panel, "Guidance circles", &showCircles);
	glui_check_circles->set_alignment(GLUI_ALIGN_RIGHT);
}
//=============================================================================
void GLViewer::initGLUIFeatures(FeatureConfig* _features, int _nFeatures)
{
	GLUI_Panel *features_panel = new GLUI_Rollout(glui, "Features", true);
	cout << _nFeatures << endl;
	features = new float[_nFeatures];
	for (int i = 0; i < _nFeatures; i++)
	{
		cout << "Feature " << i << ": " << _features[i].name << endl;
		cout << "Feature " << i << ": " << _features[i].init_value << endl;
		cout << "Feature " << i << ": " << _features[i].min_value << endl;
		cout << "Feature " << i << ": " << _features[i].max_value << endl;
		cout << "Feature " << i << ": " << _features[i].incr_value << endl;

		features[i] = _features[i].init_value;
		GLUI_Spinner *spinner =
			new GLUI_Spinner(features_panel, _features[i].name,
			&features[i], i, updateFeatures);
		spinner->set_float_limits(_features[i].min_value,
			_features[i].max_value);
		spinner->set_float_val(_features[i].init_value);
		//spinner->set_speed(_features[i].getIncrValue());
		spinner->set_alignment(GLUI_ALIGN_RIGHT);
	}
}
//=============================================================================
void GLViewer::updateFeatures(int _idxFeature)
{
	pca.editFeature(_idxFeature, features[_idxFeature]);
	pca.updateMesh(mesh);
}