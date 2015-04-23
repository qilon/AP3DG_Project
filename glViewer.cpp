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
const MyMesh::Color GLViewer::MODEL_COLOR(.9f, .79f, .69f, 1.f);

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

/* MODE BUTTON TEXT */
const char* GLViewer::GENERATE_MODE_TEXT = "Go to Generation mode";
const char* GLViewer::RECONSTRUCT_MODE_TEXT = "Go to Reconstruction mode";

/* MESH RECONSTRUCTION */
const MyMesh::Color GLViewer::SELECTED_INDEX_COLOR(1.f, .3f, .3f, 1.f);
const MyMesh::Color GLViewer::RECONSTRUCTED_POINT_COLOR(.5f, .6f, .9f, 1.f);
const MyMesh::Point GLViewer::REMOVED_POINT(0.f, 0.f, 0.f);
const int GLViewer::REMOVE_VERTEX_INDEX = 0;
const int GLViewer::REMOVE_N_RINGS = 10;
const int GLViewer::REMOVE_MAX_RINGS = 20;

//=============================================================================
/**** VARIABLES ****/

/* MAIN VARIABLES */
PCA GLViewer::pca = PCA();
MyMesh GLViewer::mesh;
float* GLViewer::features;
float* GLViewer::alphas;
int GLViewer::mode = RECONSTRUCT_MODE;
MyMesh GLViewer::recons_mesh;
VectorXi GLViewer::points_state;
int GLViewer::remove_vertex_index = REMOVE_VERTEX_INDEX;
int GLViewer::remove_n_rings = REMOVE_N_RINGS;

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
GLUI_Spinner** GLViewer::spinnersAlphas;
GLUI_Spinner** GLViewer::spinnersFeatures;
GLUI_Rollout* GLViewer::features_rollout;
GLUI_Rollout* GLViewer::alphas_rollout;
GLUI_Rollout* GLViewer::recons_rollout;
GLUI_Button* GLViewer::glui_modeButton;

//=============================================================================
void GLViewer::initGLUT(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE | GLUT_ALPHA);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - WINDOW_WIDTH) / 2,
		10);

	window_id = glutCreateWindow(WINDOW_TITLE);

	// Uncomment to enable transparencies
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_BLEND);

	glutDisplayFunc(display);
	//glutMotionFunc(motion);

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

	/* Reconstruction Panel */
	initGLUIReconstruction();

	/* Control Panel */
	initGLUIControlPanel();
}
//=============================================================================
/*
 * Function that creates a spinner for each of the features defined in pca
 */
void GLViewer::initGLUIFeatures(FeatureConfig* _features, int _nFeatures)
{
	features_rollout = new GLUI_Rollout(glui, "Features", false);
	features = new float[_nFeatures];
	spinnersFeatures = new GLUI_Spinner*[_nFeatures];
	for (int i = 0; i < _nFeatures; i++)
	{
		features[i] = _features[i].init_value;
		GLUI_Spinner *spinner =
			new GLUI_Spinner(features_rollout, _features[i].name,
			&features[i], i, updateFeature);
		spinner->set_float_limits(_features[i].min_value,
			_features[i].max_value);
		spinner->set_float_val(_features[i].init_value);
		//spinner->set_speed(_features[i].getIncrValue());
		spinner->set_alignment(GLUI_ALIGN_RIGHT);
		spinnersFeatures[i] = spinner;
	}

	if (mode == GENERATE_MODE)
	{
		features_rollout->open();
		features_rollout->enable();
	}
	else // RECONSTRUCT_MODE
	{
		features_rollout->disable();
		features_rollout->close();
	}
}
//=============================================================================
void GLViewer::initGLUIAlphas(VectorXf _alphas)
{
	alphas_rollout = new GLUI_Rollout(glui, "Alphas", false);
	alphas = new float[NUM_ALPHAS_DISPLAYED];
	spinnersAlphas = new GLUI_Spinner*[NUM_ALPHAS_DISPLAYED];
	for (int i = 0; i < NUM_ALPHAS_DISPLAYED; i++)
	{
		alphas[i] = _alphas(i);
		stringstream ss;
		ss << "Alpha No." << i << "\0";
		const string name = ss.str();
		spinnersAlphas[i] =
			new GLUI_Spinner(alphas_rollout, name.c_str(),
			&alphas[i], i, updateAlpha);
		spinnersAlphas[i]->set_float_limits(ALPHA_MIN,
			ALPHA_MAX);
		spinnersAlphas[i]->set_float_val(_alphas(i));
		//spinner->set_speed(_features[i].getIncrValue());
		spinnersAlphas[i]->set_alignment(GLUI_ALIGN_RIGHT);
	}

	if (mode == GENERATE_MODE)
	{
		alphas_rollout->open();
		alphas_rollout->enable();
	}
	else // RECONSTRUCT_MODE
	{
		alphas_rollout->disable();
		alphas_rollout->close();
	}
}
//=============================================================================
void GLViewer::initGLUIReconstruction()
{
	recons_rollout = new GLUI_Rollout(glui, "Reconstruction", false);

	GLUI_Panel* remove_panel = new GLUI_Panel(recons_rollout,"",GLUI_PANEL_RAISED);

	GLUI_Spinner* vertex_index_spinner =
		new GLUI_Spinner(remove_panel, "Vertex index",
		&remove_vertex_index, -1);
	vertex_index_spinner->set_int_val(REMOVE_VERTEX_INDEX);
	vertex_index_spinner->set_int_limits(0, recons_mesh.n_vertices()-1);
	vertex_index_spinner->set_speed(0.1f);

	GLUI_Spinner* n_rings_spinner =
		new GLUI_Spinner(remove_panel, "No. rings",
		&remove_n_rings, -1);
	n_rings_spinner->set_int_val(REMOVE_N_RINGS);
	n_rings_spinner->set_int_limits(0, REMOVE_MAX_RINGS);
	n_rings_spinner->set_speed(0.1f);

	/* Remove button */
	GLUI_Button* glui_removeButton = new GLUI_Button(remove_panel,
		"Remove points", -1, &removePointsButtonCallback);

	/* Reconstruct button */
	GLUI_Button* glui_reconstructButton = new GLUI_Button(recons_rollout,
		"Reconstruct model", -1, &reconstructButtonCallback);

	/* Reset colors button */
	GLUI_Button* glui_resetColorsButton = new GLUI_Button(recons_rollout,
		"Reset colors", -1, &resetColorsButtonCallback);

	if (mode == GENERATE_MODE)
	{
		recons_rollout->disable();
		recons_rollout->close();
	}
	else // RECONSTRUCT_MODE
	{
		recons_rollout->open();
		recons_rollout->enable();
	}
}
//=============================================================================
void GLViewer::initGLUIControlPanel()
{
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

	/* Mode button */
	glui_modeButton = glui->add_button(RECONSTRUCT_MODE_TEXT, -1, &modeButtonCallback);

	if (mode == RECONSTRUCT_MODE)
	{
		glui_modeButton->set_name(GENERATE_MODE_TEXT);
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

	// Guidance Circles
	if (showCircles)
	{
		glDisable(GL_LIGHTING);
		drawCircle(radius, 0, CIRCLE_NUM_LINES, CIRCLE_XY_COLOR);
		drawCircle(radius, 1, CIRCLE_NUM_LINES, CIRCLE_XZ_COLOR);
		drawCircle(radius, 2, CIRCLE_NUM_LINES, CIRCLE_YZ_COLOR);
		glEnable(GL_LIGHTING);
	}

	drawModel();

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
void GLViewer::modeButtonCallback(int state)
{
	updateMode();
}
//=============================================================================
void GLViewer::removePointsButtonCallback(int state)
{
	removeReconsMeshRegion(remove_vertex_index, remove_n_rings);
}
//=============================================================================
void GLViewer::reconstructButtonCallback(int state)
{
	reconstruct();
}
//=============================================================================
void GLViewer::resetColorsButtonCallback(int state)
{
	setMeshColor(recons_mesh, MODEL_COLOR);
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
	MyMesh* mesh2draw;
	glBegin(GL_TRIANGLES);
	if (mode==RECONSTRUCT_MODE)
	{
		mesh2draw = &recons_mesh;

		for (MyMesh::FaceIter f_it = mesh2draw->faces_begin();
			f_it != mesh2draw->faces_end(); ++f_it)
		{
			bool draw_triangle = true;
			MyMesh::FaceVertexIter fv_it;
			for (fv_it = mesh2draw->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				int v_idx = fv_it.handle().idx();
				draw_triangle = draw_triangle && points_state(v_idx);
			}

			if (draw_triangle)
			{
				for (fv_it = mesh2draw->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
				{
					MyMesh::Point p = mesh2draw->point(*fv_it);
					float point[3] {p[0], p[1], p[2]};

					MyMesh::Normal n = mesh2draw->normal(*fv_it);
					float normal[3] {n[0], n[1], n[2]};

					int v_idx = fv_it.handle().idx();
					MyMesh::Color c;
					if (v_idx==remove_vertex_index)
					{
						c = SELECTED_INDEX_COLOR;
					}
					else
					{
						c = mesh2draw->color(*fv_it);
					}
					glColor4f(c[0], c[1], c[2], c[3]);

					glNormal3fv(normal);
					glVertex3fv(point);
				}
			}
		}
	}
	else // GENERATE_MODE
	{
		mesh2draw = &mesh;

		for (MyMesh::FaceIter f_it = mesh2draw->faces_begin();
			f_it != mesh2draw->faces_end(); ++f_it)
		{
			MyMesh::FaceVertexIter fv_it;
			for (fv_it = mesh2draw->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				MyMesh::Point p = mesh2draw->point(*fv_it);
				float point[3] {p[0], p[1], p[2]};

				MyMesh::Normal n = mesh2draw->normal(*fv_it);
				float normal[3] {n[0], n[1], n[2]};

				MyMesh::Color c = mesh2draw->color(*fv_it);
				glColor4f(c[0], c[1], c[2], c[3]);

				glNormal3fv(normal);
				glVertex3fv(point);
			}
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
	VectorXf alphasVector = pca.getAlphas();
	for (int i = 0; i < NUM_ALPHAS_DISPLAYED; i++) {
		alphas[i] = alphasVector(i);
		spinnersAlphas[i]->set_float_val(alphasVector(i));
	}
	pca.updateMesh(mesh);
}
//=============================================================================
void GLViewer::updateAlpha(int _idxAlpha)
{
	pca.editAlpha(_idxAlpha, alphas[_idxAlpha]);
	VectorXf featuresVector = pca.getFeatures();
	for (int i = 0; i < featuresVector.size() - 1; i++) {
		features[i] = featuresVector(i);
		spinnersFeatures[i]->set_float_val(featuresVector(i));
	}
	pca.updateMesh(mesh);
}
//=============================================================================
void GLViewer::updateMode()
{
	if (mode==GENERATE_MODE)
	{
		mode = RECONSTRUCT_MODE;
		features_rollout->disable();
		features_rollout->close();
		alphas_rollout->disable();
		alphas_rollout->close();
		recons_rollout->open();
		recons_rollout->enable();
		glui_modeButton->set_name(GENERATE_MODE_TEXT);
		glui_modeButton->update_size();
	}
	else
	{
		mode = GENERATE_MODE;
		features_rollout->open();
		features_rollout->enable();
		alphas_rollout->open();
		alphas_rollout->enable();
		recons_rollout->disable();
		recons_rollout->close();
		glui_modeButton->set_name(RECONSTRUCT_MODE_TEXT);
		glui_modeButton->update_size();
	}
}
//=============================================================================
void GLViewer::removeReconsMeshRegion(int _vertex_idx, int _n_rings)
{
	// FIFO queue containing indexes of vertices to mark as unknown
	queue<int> q_vertex_idx;
	q_vertex_idx.push(_vertex_idx);

	// Mark origin vertex as unknown and change its color
	recons_mesh.point(recons_mesh.vertex_handle(_vertex_idx)) = REMOVED_POINT;
	points_state(_vertex_idx) = UNKNOWN_POINT;

	// Set of points already visited
	set<int> s_visited_idx;
	s_visited_idx.insert(_vertex_idx);

	// Mark as unknowns all the points in the _n_rings of the origin point
	for (int i = 0; i < _n_rings; i++)
	{
		int n_points_prev_ring = q_vertex_idx.size();

		for (int j = 0; j < n_points_prev_ring; j++)
		{
			int v_idx = q_vertex_idx.front();
			q_vertex_idx.pop(); // Remove from queue

			MyMesh::VertexHandle v_handle = recons_mesh.vertex_handle(v_idx);
			MyMesh::VertexVertexIter vv_it;
			for (vv_it = mesh.vv_iter(v_handle); vv_it; ++vv_it)
			{
				int v_idx = vv_it.handle().idx();
				if (s_visited_idx.count(v_idx)==0)
				{
					recons_mesh.point(*vv_it) = REMOVED_POINT;
					points_state(v_idx) = UNKNOWN_POINT;
					s_visited_idx.insert(v_idx);
					q_vertex_idx.push(v_idx);
				}
			}
		}
	}

	glutPostRedisplay();
}
//=============================================================================
void GLViewer::reconstruct()
{
	/* Reconstruct model */
	pca.reconstructMesh(recons_mesh, points_state, RECONSTRUCTED_POINT_COLOR);

	/* Set all points to known */
	points_state = VectorXi::Ones(recons_mesh.n_vertices());
}
//=============================================================================
void GLViewer::setMeshColor(MyMesh& _mesh, const MyMesh::Color& _color)
{
	MyMesh::VertexIter v_it, v_end(_mesh.vertices_end());
	for (v_it = _mesh.vertices_begin(); v_it != v_end; v_it++)
	{
		_mesh.set_color(*v_it, _color);
	}
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

	// Add vertex colors
	mesh.request_vertex_colors();

	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(mesh.vertices_end());
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		mesh.set_color(*v_it, MODEL_COLOR);
	}

	// Copy mesh to mesh to reconstruct
	recons_mesh = mesh;

	// All the points are known points
	points_state = VectorXi::Ones(recons_mesh.n_vertices())*KNOWN_POINT;

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
void GLViewer::destroy()
{
	delete[] spinnersAlphas;
	delete[] spinnersFeatures;
	delete[] alphas;
	delete[] features;
}
//=============================================================================
