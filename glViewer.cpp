#include "GLViewer.h"
//=============================================================================
int GLViewer::moving = 0;
int GLViewer::beginx = 0;
int GLViewer::beginy = 0;

GLfloat GLViewer::eye[3] = { 0.0f, 0.0f, 5.0f };
GLfloat GLViewer::center[3] = { 0.0f, 0.0f, 0.0f };
GLfloat GLViewer::fieldView = 40.0;
GLfloat GLViewer::aspectRatio = 1.0f;
GLfloat GLViewer::depthNear = 1.0f;
GLfloat GLViewer::depthFar = 50.0f;

GLfloat GLViewer::anglex = 0;   /* in degrees */
GLfloat GLViewer::angley = 0;   /* in degrees */

GLfloat GLViewer::light_diffuse[4] = { 1.0, 0.89, 0.79, 1.0 };  /* Red diffuse light. */
GLfloat GLViewer::light_ambient[4] = { 0.1, 0.1, 0.1, 1.0 };  /* Grey ambient light. */
GLfloat GLViewer::light_position[4] = { 1.0, 1.0, 1.0, 0.0 };  /* Infinite light location. */

GLfloat GLViewer::background_colour[4] = { 0.4f, 0.4f, 0.4f, 0.0f };	/* Background colour. */

PolyMesh GLViewer::mesh;

//=============================================================================
void GLViewer::initialize(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow(WINDOW_TITLE);

	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	glutKeyboardFunc(Key);
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

	gluLookAt(eye[0], eye[1], eye[2],  /* eye is at (0,0,5) */
		center[0], center[1], center[2],      /* center is at (0,0,0) */
		0.0, 1.0, 0.0);      /* up is in positive Y direction */

	/* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glClearColor(background_colour[0], background_colour[1], background_colour[2],
		background_colour[3]);

	/* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);
}
//=============================================================================
void GLViewer::display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* Adjust cube position according to user input */
	/* Note, this is very basic, and suffers from Gimball lock */
	glLoadIdentity();

	gluLookAt(eye[0], eye[1], eye[2],  /* eye is at (0,0,5) */
		center[0], center[1], center[2], /* center is at (0,0,0) */
		0.0, 1.0, 0.);      /* up is in positive Y direction */

	//glTranslatef(0.0, 0.0, -3.5);
	glRotatef(angley, 1.0, 0.0, 0.0);
	glRotatef(-anglex, 0.0, 0.0, 1.0);

	glMatrixMode(GL_MODELVIEW);

	/* Coursework: replace the simple robot. */
	drawModel();

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
	for (PolyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		PolyMesh::FaceVertexIter fv_it;
		for (fv_it = mesh.fv_iter(*f_it); fv_it; ++fv_it)
		{
			PolyMesh::Point p = mesh.point(*fv_it);
			float point[3] {p[0], p[1], p[2]};

			PolyMesh::Normal n = mesh.normal(*fv_it);
			float normal[3] {n[0], n[1], n[2]};

			glNormal3fv(normal);
			glVertex3fv(point);
		}
	}
	glEnd();
}
//=============================================================================
void GLViewer::Key(unsigned char key, int x, int y) {
	GLfloat MOVEMENT = 0.1f;
	GLfloat ANGLE = M_PI / 36;
	GLfloat COS_ANGLE = cos(ANGLE);
	GLfloat SIN_ANGLE = sin(ANGLE);
	GLfloat eye_x = eye[0];
	GLfloat eye_z = eye[2];
	switch (key) {
	case UPPER_W:
		eye[1] -= MOVEMENT;
		center[1] -= MOVEMENT;
		break;
	case LOWER_W:
		eye[1] -= MOVEMENT;
		center[1] -= MOVEMENT;
		break;
	case UPPER_S:
		eye[1] += MOVEMENT;
		center[1] += MOVEMENT;
		break;
	case LOWER_S:
		eye[1] += MOVEMENT;
		center[1] += MOVEMENT;
		break;
	case UPPER_A:
		eye[0] += MOVEMENT;
		center[0] += MOVEMENT;
		break;
	case LOWER_A:
		eye[0] += MOVEMENT;
		center[0] += MOVEMENT;
		break;
	case UPPER_D:
		eye[0] -= MOVEMENT;
		center[0] -= MOVEMENT;
		break;
	case LOWER_D:
		eye[0] -= MOVEMENT;
		center[0] -= MOVEMENT;
		break;
	case UPPER_Q:
		eye[0] = eye_x*COS_ANGLE + eye_z*SIN_ANGLE;
		eye[2] = -eye_x*SIN_ANGLE + eye_z*COS_ANGLE;
		break;
	case LOWER_Q:
		eye[0] = eye_x*COS_ANGLE + eye_z*SIN_ANGLE;
		eye[2] = -eye_x*SIN_ANGLE + eye_z*COS_ANGLE;
		break;
	case UPPER_E:
		eye[0] = eye_x*COS_ANGLE - eye_z*SIN_ANGLE;
		eye[2] = eye_x*SIN_ANGLE + eye_z*COS_ANGLE;
		break;
	case LOWER_E:
		eye[0] = eye_x*COS_ANGLE - eye_z*SIN_ANGLE;
		eye[2] = eye_x*SIN_ANGLE + eye_z*COS_ANGLE;
		break;
	case UPPER_R:
		eye[2] -= MOVEMENT;
		break;
	case LOWER_R:
		eye[2] -= MOVEMENT;
		break;
	case UPPER_F:
		eye[2] += MOVEMENT;
		break;
	case LOWER_F:
		eye[2] += MOVEMENT;
		break;
	default:
		break;
	}
	glutPostRedisplay();
}
//=============================================================================
void GLViewer::setMesh(PolyMesh& _mesh)
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
}
//=============================================================================
void GLViewer::run()
{
	glutMainLoop();
}
//=============================================================================