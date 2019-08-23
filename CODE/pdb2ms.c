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
 *   pdb2ms
 *      Converts pdb format to proper format for xms program
 *      Atom type index (defining the radius used from ms.rad file) is added to the records
 *
 *   Input:
 *      argv[1] - dry protein
 *   Output:
 *      argv[2] (or stdout) - ms file of the hydrogenless protein
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"
#define MAX_RADIUS 1000

extern void readPDB10(char *, int *, PDB **);
extern void initVDW(char *, int, PDB *);

void main(int argc, char *argv[])
{
  int	numAtoms;
  PDB	*pdbP;
  FILE *fp;
  int	it;
  int	i;
  REAL  radius[MAX_RADIUS];
  int num_radii=0;
  REAL aaa;
  char aname[256];

  if (argc < 3 || argc > 4)  {
      fprintf(stderr,"USAGE: pdb2ms protein.pdb output.ms [working directory]\n");
      exit(1);
  }

  readPDB10(argv[1], &numAtoms, &pdbP);
  /* read nonbonded parameters */
  initVDW (argv[1], numAtoms, pdbP);

  /* create table of atomic radii */
  for( i=0; i<numAtoms; i++) {
      aaa = pdbP[i].ms_rad;
      if (aaa == 0.) continue;
      for (it=0;it<num_radii;it++) { if (aaa==radius[it]) break; }
      if (it==num_radii) {
	  radius[it]=aaa; num_radii++;
      }
  }
  /* write out the table of radii */
  if (argc == 4) strcpy(aname,argv[3]);
  else strcpy(aname,".");

  strcat(aname,"/ms.rad");
  fp = fopen (aname,"w");
  if (!fp)  {
      fprintf(stderr,"*** Unable to open %s for writing\n",aname);
      exit(1);
  }
  for (i=0;i<num_radii;i++) {
      fprintf (fp,"%5d %9.5f\n",i+1,radius[i]);
  }
  fclose (fp);
/*
123451234567890
    1   2.89000
*/

  /* write out the atoms */
  fp = fopen(argv[2],"w");
  if (!fp)  {
      fprintf(stderr,"*** Unable to open %s for writing\n",argv[2]);
      exit(1);
  }

  for( i=0; i<numAtoms; i++) {
      aaa = pdbP[i].ms_rad;
      if (aaa == 0.) continue;
      for (it=0;it<num_radii;it++) { if (aaa==radius[it]) break; }
      fprintf (fp,"%10.5f%10.5f%10.5f%5d%5d%5d %s %s%s\n", 
	pdbP[i].XX, pdbP[i].YY, pdbP[i].ZZ,
	it+1,2,0,pdbP[i].resName,pdbP[i].atomType,
        pdbP[i].atomLoc);
  }
}  /* end main */
