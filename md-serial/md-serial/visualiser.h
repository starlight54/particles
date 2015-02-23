#ifndef VISUALISER_H
#define VISUALISER_H

#include <GL/freeglut.h>
#include <math.h>
#include "particleSystem.h"

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
	ParticleSystem* particles;
	double* vel;
	static void Display();
	void Init();
	void Draw();
	void GetQuaternion(double Ax, double Ay, double Az, double Bx, 
		double By, double Bz);
};

#endif