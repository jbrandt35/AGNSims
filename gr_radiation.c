

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

void rebx_gr_radiation(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){


	double* cp = rebx_get_param(sim->extras, force->ap, "c");
	if (cp == NULL){
	rebx_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
	}
	const double c = (*cp);
	struct reb_particle* const particlesl = sim->particles;

	int* ppart1 = rebx_get_param(sim->extras, force->ap, "gr_rad_part1");
	int* ppart2 = rebx_get_param(sim->extras, force->ap, "gr_rad_part2");

	if (ppart1 == NULL || ppart2==NULL){
	rebx_error(sim, "Need to set which particles to use. Set parameters: \"gr_rad_part1\" and \"gr_rad_part2\".\n");
	}

	int part1=*ppart1;
	int part2=*ppart2;
	double m0,m1;
	double rx,ry,rz;
	double rr;
	double rurx,rury,rurz;
	double rvx,rvy,rvz;
	double rv;
	double dpurv;
	double a0x,a0y,a0z;
	double a1x,a1y,a1z;
    double a0x1,a0y1,a0z1;
    double a1x1,a1y1,a1z1;
	double pfact0,pfact1;
    double n12v1,n12v2,v12,v22,v1v2;

	m0=particlesl[part1].m;
	m1=particlesl[part2].m;

	rx=particlesl[part2].x-particlesl[part1].x;
	ry=particlesl[part2].y-particlesl[part1].y;
	rz=particlesl[part2].z-particlesl[part1].z;
	rr=sqrt((rx*rx)+(ry*ry)+(rz*rz));
	rurx=rx/rr;
	rury=ry/rr;
	rurz=rz/rr;

	rvx=particlesl[part2].vx-particlesl[part1].vx;
	rvy=particlesl[part2].vy-particlesl[part1].vy;
	rvz=particlesl[part2].vz-particlesl[part1].vz;
	rv=sqrt((rvx*rvx)+(rvy*rvy)+(rvz*rvz));

	dpurv=(rvx*rurx)+(rvy*rury)+(rvz*rurz);

	double R_g=(sim->G)*m0/c;
	if (rr < 100*R_g) {

    // the following is PH2.5 (GW)
	pfact0=4*(sim->G)*(sim->G)*m0*m1*(m1/(m0+m1))/(5*pow(c,5)*pow(rr,3));
	a0x=pfact0;a0y=pfact0;a0z=pfact0;
	a0x*=((rurx*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))+(rvx*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));
	a0y*=((rury*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))+(rvy*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));
	a0z*=((rurz*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))+(rvz*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));

	pfact1=4*(sim->G)*(sim->G)*m0*m1*(m0/(m0+m1))/(5*pow(c,5)*pow(rr,3));
	a1x=pfact1;a1y=pfact1;a1z=pfact1;
	a1x*=((-rurx*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))-(rvx*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));
	a1y*=((-rury*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))-(rvy*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));
	a1z*=((-rurz*(dpurv)*((34.0*(sim->G)*(m0+m1)/(3*rr))+(6*rv*rv)))-(rvz*((-6*(sim->G)*(m0+m1)/rr)-(2*rv*rv))));

	particles[part1].ax += -a0x;particles[part1].ay += -a0y;particles[part1].az += -a0z;
    particles[part2].ax += -a1x;particles[part2].ay += -a1y;particles[part2].az += -a1z;

    // the following is PH1 (orbital precession)
    n12v1=-rurx*particlesl[part1].vx-rury*particlesl[part1].vy-rurz*particlesl[part1].vz;
    n12v2=-rurx*particlesl[part2].vx-rury*particlesl[part2].vy-rurz*particlesl[part2].vz;

    v12=particlesl[part1].vx*particlesl[part1].vx+particlesl[part1].vy*particlesl[part1].vy+particlesl[part1].vz*particlesl[part1].vz;
    v22=particlesl[part2].vx*particlesl[part2].vx+particlesl[part2].vy*particlesl[part2].vy+particlesl[part2].vz*particlesl[part2].vz;
    v1v2=particlesl[part1].vx*particlesl[part2].vx+particlesl[part1].vy*particlesl[part2].vy+particlesl[part1].vz*particlesl[part2].vz;


    a0x1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m1*m1/pow(rr,3) + (sim->G)*m1/rr/rr*(3/2*(n12v2)*(n12v2)-v12+4*v1v2-2*v22))*(-rurx) + (sim->G)*m1/rr/rr*(-4*(n12v1)+3*n12v2)*(-rvx));
    a0y1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m1*m1/pow(rr,3) + (sim->G)*m1/rr/rr*(3/2*(n12v2)*(n12v2)-v12+4*v1v2-2*v22))*(-rury) + (sim->G)*m1/rr/rr*(-4*(n12v1)+3*n12v2)*(-rvx));
    a0z1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m1*m1/pow(rr,3) + (sim->G)*m1/rr/rr*(3/2*(n12v2)*(n12v2)-v12+4*v1v2-2*v22))*(-rurz) + (sim->G)*m1/rr/rr*(-4*(n12v1)+3*n12v2)*(-rvx));

    a1x1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m0*m0/pow(rr,3) + (sim->G)*m0/rr/rr*(3/2*(n12v1)*(n12v1)-v22+4*v1v2-2*v12))*(rurx) + (sim->G)*m0/rr/rr*(4*(n12v2)-3*n12v1)*(rvx));
    a1y1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m0*m0/pow(rr,3) + (sim->G)*m0/rr/rr*(3/2*(n12v1)*(n12v1)-v22+4*v1v2-2*v12))*(rury) + (sim->G)*m0/rr/rr*(4*(n12v2)-3*n12v1)*(rvy));
    a1z1=1/c/c*((5*(sim->G)*(sim->G)*m0*m1/pow(rr,3) + 4*(sim->G)*(sim->G)*m0*m0/pow(rr,3) + (sim->G)*m0/rr/rr*(3/2*(n12v1)*(n12v1)-v22+4*v1v2-2*v12))*(rurz) + (sim->G)*m0/rr/rr*(4*(n12v2)-3*n12v1)*(rvz));
    particles[part1].ax += -a0x1;particles[part1].ay += -a0y1;particles[part1].az += -a0z1;
    particles[part2].ax += -a1x1;particles[part2].ay += -a1y1;particles[part2].az += -a1z1;

	}

}
