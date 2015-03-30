#include "glViewer.h"
//=============================================================================
void drawModel(PolyMesh& mesh)
{
	//static const int POSITION_TAG = 0;
	//static const int NORMAL_TAG   = 2;

	//glEnableVertexAttribArray(POSITION_TAG);
	//glEnableVertexAttribArray(NORMAL_TAG);

	//glVertexAttribPointer(POSITION_TAG, 3, GL_FLOAT, GL_FALSE, 0, model.vertex.data());
	//glVertexAttribPointer(NORMAL_TAG  , 3, GL_FLOAT, GL_TRUE , 0, model.normal.data());

	//typedef map< string, vector<unsigned short> >::const_iterator ConstFaceIterator;
	//for (ConstFaceIterator it = model.faces.begin(); it != model.faces.end(); ++it)
	//	glDrawElements(GL_TRIANGLES, it->second.size(), GL_UNSIGNED_SHORT, it->second.data());

	//glDisableVertexAttribArray(NORMAL_TAG);
	//glDisableVertexAttribArray(POSITION_TAG);

	glBegin(GL_TRIANGLES);
	for (PolyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		PolyMesh::FaceVertexIter fv_it;
		for (fv_it = mesh.fv_iter(*f_it); fv_it; ++fv_it)
		{
			PolyMesh::Point p = mesh.point(*fv_it);
			float point[3] {p[0], p[1], p[2]};

			//glNormal3fv(&model.normal[index * 3]);
			glVertex3fv(point);
		}
	}
	glEnd();
}
//=============================================================================
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* Adjust cube position according to user input */
	/* Note, this is very basic, and suffers from Gimball lock */
	glLoadIdentity();

	gluLookAt(eye[0], eye[1], eye[2],  /* eye is at (0,0,5) */
		center[0], center[1], center[2], /* center is at (0,0,0) */
		0.0, 1.0, 0.);      /* up is in positive Y direction */

	glTranslatef(0.0, 0.0, -3.5);
	glRotatef(angley, 1.0, 0.0, 0.0);
	glRotatef(-anglex, 0.0, 0.0, 1.0);

	glMatrixMode(GL_MODELVIEW);

	/* Coursework: replace the simple robot. */
	//robot->draw();
	drawModel(mesh);

	glutSwapBuffers();
}
//=============================================================================
void init(void)
{
	glewInit();

	/* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);
	gluPerspective( /* field of view in degree */ 40.0,
		/* aspect ratio */ 1.0,
		/* Z near */ 1.0, /* Z far */ 10.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(0.0, 0.0, 5.0,  /* eye is at (0,0,5) */
		0.0, 0.0, 0.0,      /* center is at (0,0,0) */
		0.0, 1.0, 0.);      /* up is in positive Y direction */

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
void animate(void)
{
	//cycle through color to demonstrate that animation is working
	light_diffuse[1] += 0.001f;
	if (light_diffuse[1] > 0.999f) light_diffuse[1] = 0.0f;
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);

	/* Coursework: you can remove the color animation and instead
	use it to animate your robot */
	//robot->waveHands();

	glutPostRedisplay();
}
//=============================================================================
void mouse(int button, int state, int x, int y)
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
void motion(int x, int y)
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
static void Key(unsigned char key, int x, int y) {
	GLfloat MOVEMENT = 0.1f;
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
		eye[2] += MOVEMENT;
		break;
	case LOWER_Q:
		eye[2] += MOVEMENT;
		break;
	case UPPER_E:
		eye[2] -= MOVEMENT;
		break;
	case LOWER_E:
		eye[2] -= MOVEMENT;
		break;
		//case UPPER_R:
		//	robot->startWaveLeftHand();
		//	break;
		//case LOWER_R:
		//	robot->startWaveLeftHand();
		//	break;
		//case UPPER_T:
		//	robot->startWaveRightHand();
		//	break;
		//case LOWER_T:
		//	robot->startWaveRightHand();
		//	break;
		//case UPPER_M:
		//	robot->switchModel();
		//	break;
		//case LOWER_M:
		//	robot->switchModel();
		//	break;
	default:
		break;
	}
}
//=============================================================================
void loadMesh()
{
	OpenMesh::IO::Options ropt, wopt;

	if (READ_VERTEX_COLORS)
	{
		mesh.request_vertex_colors();
		ropt += OpenMesh::IO::Options::VertexColor;
	}
	wopt += OpenMesh::IO::Options::VertexColor;

	cout << "Reading file... " << endl;
	if (!OpenMesh::IO::read_mesh(mesh, PLY_FILENAME))
	{
		cout << "Could not read file: " << PLY_FILENAME << endl << endl;

		cout << "Press any key to exit... ";
		getchar();

		//return -1;
	}
}
//=============================================================================
void initViewer(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutCreateWindow("AP3DG");

	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	glutKeyboardFunc(Key);
	glutMotionFunc(motion);
	glutIdleFunc(animate);

	init();
}
//=============================================================================