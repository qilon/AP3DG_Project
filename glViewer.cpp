#include "GLViewer.h"
//=============================================================================
const char* GLViewer::WINDOW_TITLE = "AP3DG";
const int GLViewer::WINDOW_WIDTH = 800;
const int GLViewer::WINDOW_HEIGHT = 800;

const GLfloat GLViewer::EYE_ZOOM_INCR = 0.1f;
const GLfloat GLViewer::EYE_ANGLE_INCR = M_PI / 36;

const GLfloat GLViewer::RADIUS_OFFSET = .1f;
const GLint GLViewer::CIRCLE_NUM_LINES = 100;
const Vector3f GLViewer::CIRCLE_XY_COLOR = {1.f, 0.f, 0.f};
const Vector3f GLViewer::CIRCLE_XZ_COLOR = { 0.f, 1.f, 0.f };
const Vector3f GLViewer::CIRCLE_YZ_COLOR = { 0.f, 0.f, 1.f };

const Vector3f GLViewer::MODEL_COLOR = { 1.f, 0.89f, 0.79f };  /* Diffuse light. */

int GLViewer::moving = 0;
int GLViewer::beginx = 0;
int GLViewer::beginy = 0;

Vector3f GLViewer::eye = { 0.0f, 0.0f, eyeDistance };
GLfloat GLViewer::eyeDistance = 5.0f;
Matrix3f GLViewer::eyeRotation = Matrix3f::Identity();
Vector3f GLViewer::center = { 0.0f, 0.0f, 0.0f };
Vector3f GLViewer::up = { 0.0f, 1.0f, 0.0f };

GLfloat GLViewer::fieldView = 40.0;
GLfloat GLViewer::aspectRatio = 1.0f;
GLfloat GLViewer::depthNear = 1.0f;
GLfloat GLViewer::depthFar = 50.0f;

GLfloat GLViewer::anglex = 0;   /* in degrees */
GLfloat GLViewer::angley = 0;   /* in degrees */

GLfloat GLViewer::light_ambient[4] = { 0.1, 0.1, 0.1, 1.0 };  /* Ambient light. */
GLfloat GLViewer::light_position[4] = { 1.0, 1.0, 1.0, 0.0 };  /* Infinite light location. */

GLfloat GLViewer::background_colour[4] = { 0.4f, 0.4f, 0.4f, 0.0f };	/* Background colour. */

MyMesh GLViewer::mesh;

GLfloat GLViewer::radius = 1.f;
bool GLViewer::showCircles = true;

int GLViewer::idxFeature = 0;

PCA GLViewer::pca = PCA();

GLUI* GLViewer::glui1;
GLUI* GLViewer::glui2;
int GLViewer::window_id;

//=============================================================================
void GLViewer::myGlutIdle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != window_id)
		glutSetWindow(window_id);

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::initialize(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	int window_id = glutCreateWindow(WINDOW_TITLE);

	glui1 = GLUI_Master.create_glui_subwindow(window_id,
		GLUI_SUBWINDOW_RIGHT);
	glui2 = GLUI_Master.create_glui_subwindow(window_id,
		GLUI_SUBWINDOW_BOTTOM);


	glui1->add_checkbox("Lighting");

	glui1->set_main_gfx_window(window_id);



	GLUI_Master.set_glutIdleFunc(myGlutIdle);

	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	//glutKeyboardFunc(Key);
	GLUI_Master.set_glutKeyboardFunc(Key);
	glutMotionFunc(motion);

	init();
}
//=============================================================================
void GLViewer::init(void)
{
	glewInit();

	/* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);
	gluPerspective( /* field of view in degree */ fieldView,
		/* aspect ratio */ aspectRatio,
		/* Z near */ depthNear, /* Z far */ depthFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(eye(0), eye(1), eye(2),	/* eye */
		center(0), center(1), center(2),	/* center */
		up(0), up(1), up(2));      /* up is in positive Y direction */

	/* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glEnable(GL_COLOR_MATERIAL);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glClearColor(background_colour[0], background_colour[1], 
		background_colour[2], background_colour[3]);

	/* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);

	/* Anti-aliasing*/
	glEnable(GL_MULTISAMPLE);
}
//=============================================================================
void GLViewer::display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* Adjust cube position according to user input */
	/* Note, this is very basic, and suffers from Gimball lock */
	glLoadIdentity();

	gluLookAt(eye[0], eye[1], eye[2],  /* eye */
		center[0], center[1], center[2], /* center */
		up[0], up[1], up[2]);      /* up is in positive Y direction */

	//glTranslatef(meshTranslation[0], meshTranslation[1], meshTranslation[2]);
	//glRotatef(anglex, 1.0, 0.0, 0.0);
	//glRotatef(angley, 0.0, 1.0, 0.0);
	//glRotatef(anglez, 0.0, 0.0, 1.0);

	glMatrixMode(GL_MODELVIEW);

	drawModel();

	if (showCircles)
	{
		drawCircle(radius, center, 0, CIRCLE_NUM_LINES, CIRCLE_XY_COLOR);
		drawCircle(radius, center, 1, CIRCLE_NUM_LINES, CIRCLE_XZ_COLOR);
		drawCircle(radius, center, 2, CIRCLE_NUM_LINES, CIRCLE_YZ_COLOR);
	}

	// Draw text with information on eigenvectors
	stringstream ss;
	glColor3f(0, 0, 0);
	ss << "Feature number: " << idxFeature;
	string str = ss.str();
	drawText(str.data(), str.size(), 0, 590);
	ss.str("");
	ss.clear();
	ss << "Value: " << pca.getFeature(idxFeature);
	str = ss.str();
	drawText(str.data(), str.size(), 0, 580);

	glutSwapBuffers();
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
}
//=============================================================================
void GLViewer::motion(int x, int y)
{
	if (moving) {
		anglex = anglex + (x - beginx);
		angley = angley + (y - beginy);
		beginx = x;
		beginy = y;
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
void GLViewer::Key(unsigned char key, int x, int y) {
	switch (key) {
	case UPPER_W:
		rotateEye(-EYE_ANGLE_INCR, 0);
		break;
	case LOWER_W:
		rotateEye(-EYE_ANGLE_INCR, 0);
		break;
	case UPPER_S:
		rotateEye(EYE_ANGLE_INCR, 0);
		break;
	case LOWER_S:
		rotateEye(EYE_ANGLE_INCR, 0);
		break;
	case UPPER_A:
		rotateEye(-EYE_ANGLE_INCR, 1);
		break;
	case LOWER_A:
		rotateEye(-EYE_ANGLE_INCR, 1);
		break;
	case UPPER_D:
		rotateEye(EYE_ANGLE_INCR, 1);
		break;
	case LOWER_D:
		rotateEye(EYE_ANGLE_INCR, 1);
		break;
	case UPPER_Q:
		rotateEye(-EYE_ANGLE_INCR, 2);
		break;
	case LOWER_Q:
		rotateEye(-EYE_ANGLE_INCR, 2);
		break;
	case UPPER_E:
		rotateEye(EYE_ANGLE_INCR, 2);
		break;
	case LOWER_E:
		rotateEye(EYE_ANGLE_INCR, 2);
		break;
	case UPPER_R:
		zoomEye(-EYE_ZOOM_INCR);
		break;
	case LOWER_R:
		zoomEye(-EYE_ZOOM_INCR);
		break;
	case UPPER_F:
		zoomEye(EYE_ZOOM_INCR);
		break;
	case LOWER_F:
		zoomEye(EYE_ZOOM_INCR);
		break;
	case UPPER_C:
		showCircles = !showCircles;
		break;
	case LOWER_C:
		showCircles = !showCircles;
		break;
	case UPPER_K:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) - 1.f);
		pca.updateMesh(mesh);
		break;
	case LOWER_K:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) - 1.f);
		pca.updateMesh(mesh);
		break;
	case UPPER_L:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) + 1.f);
		pca.updateMesh(mesh);
		break;
	case LOWER_L:
		pca.editFeature(idxFeature, pca.getFeature(idxFeature) + 1.f);
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
		if (idxFeature < pca.getControllers() - 1) {
			idxFeature++;
		}
		break;
	case LOWER_O:
		if (idxFeature < pca.getControllers() - 1) {
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

	updateCenterEye();

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::rotateEye(GLfloat angle, int axis)
{
	GLfloat cos_angle = cos(angle);
	GLfloat sin_angle = sin(angle);

	Vector3f originalEye;
	originalEye << 0.f, 0.f, eyeDistance;

	Matrix3f newRotation = Matrix3f::Identity();

	switch (axis)
	{
	case 0: // x
		newRotation(1, 1) = cos_angle;
		newRotation(1, 2) = -sin_angle;
		newRotation(2, 1) = sin_angle;
		newRotation(2, 2) = cos_angle;
		break;
	case 1: // y
		newRotation(0, 0) = cos_angle;
		newRotation(0, 2) = sin_angle;
		newRotation(2, 0) = -sin_angle;
		newRotation(2, 2) = cos_angle;
		break;
	case 2: // z
		newRotation(0, 0) = cos_angle;
		newRotation(0, 1) = -sin_angle;
		newRotation(1, 0) = sin_angle;
		newRotation(1, 1) = cos_angle;
		break;
	default:
		break;
	}

	eyeRotation = newRotation * eyeRotation;

	eye = center + (eyeRotation * originalEye);

	up = eyeRotation.col(1);
}
//=============================================================================
void GLViewer::zoomEye(GLfloat distance)
{
	eyeDistance = eyeDistance + distance;
	Vector3f originalEye;
	originalEye << 0.f, 0.f, eyeDistance;

	eye = center + eyeRotation * originalEye;
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
	// Translation
	glTranslatef(center(0), center(1), center(2));

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
	pca.read(_pca_filename_url, _features_filename_url);

	// Update mesh based on the pca
	pca.updateMesh(mesh);

	updateCenterEye();
}
//=============================================================================
void GLViewer::updateCenterEye()
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

	center(0) = min_x + (max_x - min_x) / 2;
	center(1) = min_y + (max_y - min_y) / 2;
	center(2) = min_z + (max_z - min_z) / 2;

	eye(0) = center(0);
	eye(1) = center(1);
	eye(2) = center(2) + eyeDistance;

	up(0) = 0.0f;
	up(1) = 1.0f;
	up(2) = 0.0f;

	radius = center.cwiseAbs().maxCoeff() + RADIUS_OFFSET;
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
