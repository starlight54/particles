#include "iSimulationIterator.h"

ISimulationIterator::ISimulationIterator()
{
	nIters = simulationTime / deltaT;
}

ISimulationIterator::~ISimulationIterator()
{

}
