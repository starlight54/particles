#include "visualiser.h"

Visualiser* Visualiser::Get()
{
	static Visualiser visualiser;
	return &visualiser;
}

void Visualiser::GlutInit(int argc, char* argv[])
{
	glutInit(&argc, argv);
	Init();
}

Visualiser::~Visualiser()
{
	gluDeleteQuadric(obj);
}

void Visualiser::SetParticles(ParticleSystem* particles)
{
	this->particles = particles;
}

void Visualiser::Init()
{
	//Create window with settings
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("MolDyn");

	//Set up OpenGL params
	static const GLfloat blackColour[4] = {0.0, 0.0, 0.0, 1.0};
	static const GLfloat whiteColour[4] = {1.0, 1.0, 1.0, 1.0};
	static const GLfloat redColour[4] = {1.0, 0.3, 0.3, 1.0};
	static const GLfloat light1_pos[] = {-100.0, 100.0, 0.0, 1.0};
	static const GLfloat light2_pos[] = {100, -100, 100, 1.0};

	glLightfv(GL_LIGHT0, GL_POSITION, light1_pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, redColour);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteColour);
	glLightfv(GL_LIGHT0, GL_SPECULAR, whiteColour);

	glLightfv(GL_LIGHT1, GL_POSITION, light2_pos);
	glLightfv(GL_LIGHT1, GL_AMBIENT, whiteColour);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, redColour);
	glLightfv(GL_LIGHT1, GL_SPECULAR, redColour);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	glClearColor(0, 0, 0, 0);
	glLineWidth(1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-100.0, 100.0, -100.0, 100.0, 100.0, -100.0);
	glRotated(25, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);

	glLightfv(GL_LIGHT0, GL_POSITION, light1_pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, redColour);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteColour);
	glLightfv(GL_LIGHT0, GL_SPECULAR, whiteColour);

	glLightfv(GL_LIGHT1, GL_POSITION, light2_pos);
	glLightfv(GL_LIGHT1, GL_AMBIENT, whiteColour);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, redColour);
	glLightfv(GL_LIGHT1, GL_SPECULAR, redColour);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	glutDisplayFunc(Visualiser::Display);
	obj = gluNewQuadric();
	gluQuadricDrawStyle(obj, GLU_FILL);
}

void Visualiser::Draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glColor3d(1, 1, 1);
	glutWireCube(100);
	glColor3d(1, 0, 0);
	if (particles) {
		for (int i = 0; i < particles->numberParticles; i++) {
			glPushMatrix();
			glTranslated(particles->pos[i * 3 + 0] - 50,
				particles->pos[i * 3 + 1] - 50, 
				particles->pos[i * 3 + 2] - 50);
			gluSphere(obj, 1.0, 10, 10);
			glPopMatrix();
		}
	}
	glutSwapBuffers();
}

void Visualiser::Display()
{
	Visualiser::Get()->Draw();
}