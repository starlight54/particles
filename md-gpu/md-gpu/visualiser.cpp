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
	static const GLfloat light1_pos[] = {-50.0, 50.0, 0.0, 1.0};
	static const GLfloat light2_pos[] = {100, -100, 100, 1.0};

	glClearColor(0, 0, 0, 0);
	glLineWidth(1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-100.0, 100.0, -100.0, 100.0, 100.0, -100.0);
	glRotated(25, 1, 1, 0);
	//glRotated(90,0,1,0);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);

	glLightfv(GL_LIGHT0, GL_POSITION, light1_pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, whiteColour);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, redColour);
	glLightfv(GL_LIGHT0, GL_SPECULAR, blackColour);

	glLightfv(GL_LIGHT1, GL_POSITION, light2_pos);
	glLightfv(GL_LIGHT1, GL_AMBIENT, whiteColour);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, redColour);
	glLightfv(GL_LIGHT1, GL_SPECULAR, redColour);

	//glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);

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
	//glutSolidCylinder()
	if (particles && vel) {
		for (int i = 0; i < particles->numParticles; i++) {
			glPushMatrix();
			glTranslated(particles->pos[i * 3 + 0] - 50,
				particles->pos[i * 3 + 1] - 50,
				particles->pos[i * 3 + 2] - 50);
			glColor3d(1, 0, 0);
			gluSphere(obj, 1.0, 10, 10);
			SetZAxisDirection(vel[i * 3 + 0], vel[i * 3 + 1],
				vel[i * 3 + 2]);
			modB = sqrt(pow(vel[i * 3 + 0], 2) +
				pow(vel[i * 3 + 1], 2) + pow(vel[i * 3 + 2], 2));
			glColor3d(0, 1, 0);
			gluCylinder(obj, 0.1, 0, modB, 10, 10);
			glPopMatrix();
		}
	}
	glutSwapBuffers();
}

void Visualiser::SetZAxisDirection(double x, double y, double z)
{
	modB = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	Wx = -x / modB;
	Wy = -y / modB;
	Wz = -z / modB;

	double rotAngle = acos(Wz);
	if (rotAngle != 0) {
		qX = Wy;
		qY = -Wx;
		qZ = 0;
		modB = sqrt(pow(qX, 2) + pow(qY, 2) + pow(qZ, 2));
		glRotated(-rotAngle/M_PI*180, qX / modB, qY / modB, qZ / modB);
	}
}

Quaternion Visualiser::GetQuaternion(double Bx, double By, double Bz, double Ax,
	double Ay, double Az)
{
	modBSquared = (pow(Bx, 2) + pow(By, 2) + pow(Bz, 2));
	normUNormV = sqrt((pow(Ax, 2) + pow(Ay, 2) + pow(Az, 2))
		* modBSquared);
	realPart = normUNormV + Ax * Bx + Ay * By + Az * Bz;

	if (realPart < 1.e-6f * normUNormV) {
		/* If u and v are exactly opposite, rotate 180 degrees
			* around an arbitrary orthogonal axis. Axis normalisation
			* can happen later, when we normalise the quaternion. */
		realPart = 0.0f;
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

	double quat_mod = sqrt(pow(realPart, 2) + pow(Wx, 2) +
		pow(Wy, 2) + pow(Wz, 2));

	qR = realPart / quat_mod;
	qX = Wx / quat_mod;
	qY = Wy / quat_mod;
	qZ = Wz / quat_mod;

	return Quaternion(qR, qX, qY, qZ);
}

void Visualiser::Display()
{
	Visualiser::Get()->Draw();
}

Quaternion::Quaternion(double real, double x, double y, double z)
{
	this->real = real;
	this->x = x;
	this->y = y;
	this->z = z;
}

Quaternion::~Quaternion()
{

}

Quaternion::Quaternion()
{

}