#ifndef VISUALISER_H
#define VISUALISER_H
#define _USE_MATH_DEFINES

#include <GL/freeglut.h>
#include <math.h>
#include "particleSystem.h"

class Quaternion
{
public:
	Quaternion();
	Quaternion(double real, double x, double y, double z);
	~Quaternion();
	double real;
	double x;
	double y;
	double z;
};

class Visualiser
{
public:
	static Visualiser* Get();
	void GlutInit(int argc, char* argv[]);
	~Visualiser();
	void SetData(ParticleSystem* particles, double* vel);
private:
	//Required for callback when drawing
	GLUquadricObj* obj;
	Quaternion quat;
	ParticleSystem* particles;
	double* vel;
	static void Display();
	void Init();
	void Draw();
	void SetZAxisDirection(double x, double y, double z);
	double qR, qX, qY, qZ, Wx, Wy, Wz, quatMod, normUNormV, realPart,
		modBSquared, modB, theta, nextTheta, x1, x2, y1, y2;
	Quaternion GetQuaternion(double Ax, double Ay, double Az, double Bx, 
		double By, double Bz);
};

#endif