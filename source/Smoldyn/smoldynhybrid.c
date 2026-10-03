/* written by Martin Robinson and Steve Andrews, started 5/2013.
 This is part of the VCell interface for the Smoldyn program.
 This is only compiled during VCell compilation and is only run from VCell code.
 See documentation called SmoldynManual.pdf and SmoldynCodeDoc.pdf, and the Smoldyn
 website, which is at www.smoldyn.org.
 Copyright 2003-2016 by Steven Andrews.  This work is distributed under the terms
 of the Gnu Lesser General Public License (LGPL). */

#include "smoldyn.h"
#include "random2.h"
#include "smoldynfuncs.h"
#include <string>
using std::string;

#include <algorithm>
using namespace std;

/* Positions hold sim->dim coordinates; pad to 3 so 1D/2D models never read past them. */
static void pad3(simptr sim, const double* pos, double* p3) {
	p3[0] = p3[1] = p3[2] = 0;
	for (int d = 0; d < sim->dim && d < 3; d++) p3[d] = pos[d];
}

double evaluateVolRnxRate(simptr sim, rxnptr reaction,  double* pos){
	double p[3];
	pad3(sim, pos, p);
	return reaction->rateValueProvider->getValue(sim->time, p[0], p[1], p[2], reaction);
}

double evaluateMemRnxRate(simptr sim, rxnptr reaction, double* pos, char* panelName)
{
	double p[3];
	pad3(sim, pos, p);
	return reaction->rateValueProvider->getValue(sim->time, p[0], p[1], p[2], reaction, panelName);
}

double evaluateSurfActionRate(simptr sim, surfactionptr actdetails, MolecState ms, double* pos, char* panelName)
{
	double p[3];
	pad3(sim, pos, p);
	return actdetails->srfRateValueProvider[ms]->getValue(sim->time, p[0], p[1], p[2], actdetails, panelName);
}


int randomPosInMesh(simptr sim, double* centerPos,  double* pos)
{
	int dim = sim -> dim;
	double deltaX, deltaY, deltaZ;
	double originX, originY, originZ;
	double domainSizeX, domainSizeY, domainSizeZ;
	double delta[3];
	sim->mesh->getDeltaXYZ(delta);

	deltaX = delta[0];
	originX = sim->wlist[0]->pos; //wlist[0],[1] are x origin and x far most 
	domainSizeX = sim->wlist[1]->pos - sim->wlist[0]->pos;
	pos[0]=unirandCCD(max(originX, (centerPos[0]-0.5*deltaX)) , min((originX+domainSizeX), (centerPos[0] + 0.5*deltaX)));
	if(dim > 1)
	{
		deltaY = delta[1];
		originY = sim->wlist[2]->pos; //wlist[2],[3] are y origin and y far most
		domainSizeY = sim->wlist[3]->pos - sim->wlist[2]->pos;
		pos[1]=unirandCCD(max(originY, (centerPos[1]-0.5*deltaY)) , min((originY+domainSizeY), (centerPos[1] + 0.5*deltaY)));
	}
	if(dim > 2)
	{
		deltaZ = delta[2];
		originZ = sim->wlist[4]->pos; //wlist[4],[5] are z origin and z far most
		domainSizeZ = sim->wlist[5]->pos - sim->wlist[4]->pos;
		pos[2]=unirandCCD(max(originZ, (centerPos[2]-0.5*deltaZ)) , min((originZ+domainSizeZ), (centerPos[2] + 0.5*deltaZ)));
	}

	return 0;
}


