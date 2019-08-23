/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * This is part of the DOWSER program
 *
 * DOWSER finds buried water molecules in proteins.
 *  cf. Zhang & Hermans, Proteins: Structure, Function and Genetics 24: 433-438, 1996
 *
 * DOWSER was developed by the Computational Structural Biology Group at the 
 * University of North Carolina, Chapel Hill by Li Zhang, Xinfu Xia,Jan Hermans, 
 * and Dave Cavanaugh.  Revised 1998.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   drain.c
 *      Remove surface waters
 *
 *   Input:
 *      argv[1] - protein structure
 *      argv[2] - water molecules
 *   Output:
 *      argv[3] - surface waters
 *      argv[4] - buried waters
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define	MAXCAGE	2000
#define TRUE 1
#define FALSE 0
#define NUM 40
#define CUTOFF 100.0
#define ANGRID 3.1415926/NUM 
#define probe 0.8 
#define D2R 3.1415926/180.0  /* DegreeToRadians */
#define	theta_step 10.0
#define	phi_step theta_step/sin(theta_r);

extern void readPDB10 (char *, int *, PDB **);
extern void WritePDB (FILE*, int, PDB *);

typedef struct cagedatom
    { REAL dx; REAL dy; REAL dz; REAL rdisc; REAL r2; }
    CagedAtom;

void main(int argc, char *argv[])
{
int	inN;
PDB	*inAtom;
int	 tarN,  objN;
PDB	*tarP, *objP;
int	i, jat;
REAL	dx, dy, dz;
REAL	dist;
int	cage_num;
CagedAtom ca[2000];
int	water_escape, collision;
REAL	theta, theta_r, phi, phi_r, projection;
XYZ	arrow;
FILE	*fp, *surfWat, *buriedWat;

    if (argc < 4 || argc > 5)  {
	fprintf(stderr, "USAGE: drain protein.pdb water.pdb surface.pdb buried.pdb (or stdout)\n\n");
	exit(1);
    }

    readPDB10 (argv[2], &inN, &inAtom);
    if (EQUAL(inAtom[1].atomType, "H")) {  /* not pure Oxygen */
	objN = inN/3;
	objP = (PDB *) malloc (objN*sizeof(PDB));
	for (i=0; i<objN; i++) objP[i] = inAtom[i*3];
    }
    else {
	objN = inN;
	objP = (PDB *) malloc (objN*sizeof(PDB));
	for (i=0; i<objN; i++) { objP[i] = inAtom[i]; }
    }

    readPDB10 (argv[1], &tarN, &tarP);  /* protein framework */

    surfWat = fopen(argv[3], "w");
    if (argc == 4)  buriedWat = stdout;
    else  buriedWat = fopen(argv[4], "w");

    for (i=0; i<tarN; i++) {
	if  (EQU(tarP[i].atomType, "c")) { /* carbon block */
	    if (EQUAL(tarP[i].atomType, "c")) tarP[i].occupancy = 1.76;
	    else if (EQUAL(tarP[i].atomType,"ca")) tarP[i].occupancy = 1.87;
	    else if (EQUAL(tarP[i].atomType,"cb")) tarP[i].occupancy = 1.87;

	    else if (EQUAL(tarP[i].resName,"tyr")    /*non-charged small**/
		|| EQUAL(tarP[i].resName,"his")
		|| EQUAL(tarP[i].resName,"phe")
		|| EQUAL(tarP[i].resName,"trp")) 
		tarP[i].occupancy = 1.76;

	    else if
		 ((EQUAL(tarP[i].resName,"asp")&& EQUAL(tarP[i].atomType,"cg"))
		||(EQUAL(tarP[i].resName,"asn")&& EQUAL(tarP[i].atomType,"cg"))
		||(EQUAL(tarP[i].resName,"glu")&& EQUAL(tarP[i].atomType,"cd"))
		||(EQUAL(tarP[i].resName,"gln")&& EQUAL(tarP[i].atomType,"cd"))
		||(EQUAL(tarP[i].resName,"arg")&& EQUAL(tarP[i].atomType,"cz")))
		tarP[i].occupancy = 1.76; 

	    else tarP[i].occupancy = 1.87;
	}
	/* non-carbon block */
	else if (EQU(tarP[i].atomType, "n")) tarP[i].occupancy = 1.7;
	else if (EQU(tarP[i].atomType, "o")) tarP[i].occupancy = 1.4;
	else if (EQU(tarP[i].atomType, "s")) tarP[i].occupancy = 1.9;
	else tarP[i].occupancy = 1.5;
    }

    for  (i=0; i<objN; i++) {
	/*make a list of atoms near the water*/
	cage_num = 0;
	for  (jat = 0; jat < tarN; jat++) {
	    dx = tarP[jat].XX - objP[i].XX;
	    dy = tarP[jat].YY - objP[i].YY;
	    dz = tarP[jat].ZZ - objP[i].ZZ;
	    dist = dx*dx + dy*dy + dz*dz;
	    if (dist<CUTOFF) {
		ca[cage_num].dx = dx;
		ca[cage_num].dy = dy;
		ca[cage_num].dz = dz;
		ca[cage_num].r2 = dist;
		ca[cage_num].rdisc =
		    (tarP[jat].occupancy+probe)*(tarP[jat].occupancy+probe);
		cage_num++;
	    }
	}
	if (cage_num>MAXCAGE) {
	    fprintf(stderr, "PAIRLIST size of %5d ", cage_num);
	    fprintf(stderr,"exceeds maximum of MAXCAGE\n");
	    exit(1);
	}

	/* check if the water can get out straight */
	/* for each direction */
	water_escape = FALSE;
	for (theta=theta_step; theta<(180.0-theta_step) && !water_escape;) {
	    theta_r = theta*D2R;
	    for (phi=0; phi<360;) {  /*for each direction**/
		phi_r = phi*D2R;
		arrow.XX = sin(theta_r)*cos(phi_r);
		arrow.YY = sin(theta_r)*sin(phi_r);
		arrow.ZZ = cos(theta_r);

		collision = FALSE;
		for  (jat=0; jat<cage_num; jat++) {
		    projection = ca[jat].dx*arrow.XX
			+ ca[jat].dy*arrow.YY
			+ ca[jat].dz*arrow.ZZ;
		    if  (projection<=0.0) continue;
		    if  (ca[jat].r2-projection*projection < ca[jat].rdisc) {  
			collision = TRUE;	/*	Collision	*/
			break;
		    }
		}
		if(!collision) {
		    water_escape = TRUE; /* This water can get out unscratched.   */
		    break;
		}
		phi += phi_step;
	    }			/*	for phi	*/
	    theta += theta_step;
	} /* for theta */

	if(water_escape) {
	    fp = surfWat;
	    strcpy (objP[i].segID, "surf");
	}
	else {
	    fp = buriedWat;
	    strcpy (objP[i].segID, "bury");
	}
	WritePDB (fp, 1, &(objP[i]));
	if (objN==inN/3) {
	    WritePDB (fp, 1, &(inAtom[3*i+1]));
	    WritePDB (fp, 1, &(inAtom[3*i+2]));
	}
    } /* for each obj */
}   /* end main() "drain" * */
/******************************************************************************/
/*									      */
/*									      */
/*	Algorithm:							      */
/*	for each water							      */
/*	  looking into all the directions				      */
/*	  for each direction						      */
/*	    checking all the protein atoms				      */
/*	    for each protein atom,					      */
/*	      checking if it will be hit by the out going water		      */
/*	      YES							      */
/*		  stop checking the remaining protein atoms		      */
/*		  go check next direction				      */
/*	      NO							      */
/*		  check the next protein atom				      */
/*	    end for each protein atom					      */
/*	    if (NO_HIT) the water is not buried, check next water	      */
/*	    else go check next direction				      */
/*	  end for each direction					      */
/*	  if (water is hit in all the directions) tag the water as buried     */
/*	end for each water						      */
/******************************************************************************/
