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

void Visualiser::SetData(ParticleSystem* particles, double* vel)
{
	this->particles = particles;
	this->vel = vel;
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

void Visualiser::GetQuaternion(double Ax, double Ay, double Az, double Bx, 
	double By, double Bz)
{
		double norm_u_norm_v = sqrt((pow(Ax, 2) + pow(Ay, 2) + pow(Az,2)) 
			* (pow(Bx, 2) + pow(By, 2) + pow(Bz, 2)));
		double real_part = norm_u_norm_v + Ax * Bx + Ay * By + Az * Bz;
		double Wx, Wy, Wz;

		if (real_part < 1.e-6f * norm_u_norm_v) {
			/* If u and v are exactly opposite, rotate 180 degrees
			* around an arbitrary orthogonal axis. Axis normalisation
			* can happen later, when we normalise the quaternion. */
			real_part = 0.0f;
			if (abs(Ax) > abs(Az)) {
				Wx = -Ay;
				Wy = Ax;
				Wz = 0.f;
			} else {
				Wx = 0.f;
				Wy = -Az;
				Wz = Ay;
			}
		} else {
			/* Otherwise, build quaternion the standard way. */
			Wx = Ay * Bz - Az * By;
			Wy = Az * Bx - Ax * Bz;
			Wz = Ax * By - Ay * Bx;
		}

		double quat_mod = sqrt(pow(real_part, 2) + pow(Wx, 2) +
			pow(Wy, 2) + pow(Wz, 2));

		double qR, qX, qY, qZ;

	
}

void Visualiser::Display()
{
	Visualiser::Get()->Draw();
}