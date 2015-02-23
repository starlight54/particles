#ifndef VISUALISER_H
#define VISUALISER_H

#include <GL/freeglut.h>
#include "particleSystem.h"

class Visualiser
{
public:
	static Visualiser* Get();
	void GlutInit(int argc, char* argv[]);
	~Visualiser();
	void SetParticles(ParticleSystem* particles);
private:
	//Required for callback when drawing
	GLUquadricObj* obj;
	ParticleSystem* particles;
	static void Display();
	void Init();
	void Draw();
};

#endif