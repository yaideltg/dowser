/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * This is part of the DOWSER program
 *
 * DOWSER finds buried water molecules in proteins.
 *  cf. Zhang & Hermans, Proteins: Structure, Function and Genetics 24: 433-438, 1996
 *
 * DOWSER was developed by the Computational Structural Biology Group at the 
 * University of North Carolina, Chapel Hill by Li Zhang, Xinfu Xia, Jan Hermans, 
 * and Dave Cavanaugh.  Revised 1998.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   readPDB10 
 *      Function to parse a pdb file
 *
 *   Input:
 *      char *pdb_file - pdb filename
 *   Output:
 *      int *numAtoms - number of atoms
 *      PDB **Atoms - atomic coordinates in a structure of type PDB
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define	MAXPDB	200000    /* Maximium number of atoms */
#define MAXCHARLINE 255

void	fld2s(char *, char * );
void	fld2i(char *, int, int * );
void	readCol(char *, int, int, char * );

void readPDB10(char *pdb_file, int *numAtoms, PDB **Atoms)
{
FILE    *file_ref;
char	line[101];
char	field00[7], field01[6], field02[2], field03[3], field04[3],
	field05[2], field06[4], field07[2], field08[2], field09[5],
	field10[2], field11[4], field12[9], field13[9], field14[9],
	field15[7], field16[7], field17[2], field18[4], field19[3],
	field20[5], field21[3], field22[3];
char	recdName[7];	/*	 1 -  6	*/

int	i = 0;
int	j, lineWidth, nAtom;
PDB	*atoms, *pureAtoms;

atoms = (PDB *) malloc ((MAXPDB)*sizeof(PDB));
if (!atoms)  {
    fprintf(stderr,"*** Unable to allocate necessary memory in readPDB10()\n");
    exit(1);
}

/* Read pdb structure from stdin or a file */
if (pdb_file == NULL)  { 
    file_ref = stdin;
}
else  {
    file_ref = fopen(pdb_file, "r");
    if (file_ref == NULL)  { 
	fprintf(stderr, "Can't open %s file!\n", pdb_file); 
	exit(1);
    }
}

while (TRUE)  {
    fgets(line, 100, file_ref);
    if (feof(file_ref))  break;
    sscanf(line, "%6s", recdName);
    if (EQUAL(recdName,"atom") || EQUAL(recdName,"hetatm"))  {
	lineWidth = 80;
	for( j=0; j<80; j++)  {
	    if (line[j]=='\n')  { lineWidth = j; break; }
	} 
	for (j=lineWidth;j<80;j++)  line[j]=' ';

	readCol(line,  1,  6, field00);  fld2s(field00, atoms[i].recdName);
	readCol(line,  7, 11, field01);  fld2i(field01, 5,&(atoms[i].serial));
	readCol(line, 12, 12, field02);
	readCol(line, 13, 14, field03);  fld2s(field03, atoms[i].atomType);
	readCol(line, 15, 16, field04);  fld2s(field04, atoms[i].atomLoc);
	readCol(line, 17, 17, field05);  fld2s(field05, atoms[i].altLoc);
	readCol(line, 18, 20, field06);  fld2s(field06, atoms[i].resName);
	readCol(line, 21, 21, field07);
	readCol(line, 22, 22, field08);  fld2s(field08, atoms[i].chainID);
	readCol(line, 23, 27, field09);  fld2i(field09, 5,&(atoms[i].resSeq));
	readCol(line, 27, 27, field10);  fld2s(field10, atoms[i].iCode);
	readCol(line, 28, 30, field11);
	readCol(line, 31, 38, field12);  atoms[i].XX = (REAL) atof(field12);
	readCol(line, 39, 46, field13);  atoms[i].YY = (REAL) atof(field13);
	readCol(line, 47, 54, field14);  atoms[i].ZZ = (REAL) atof(field14);
	readCol(line, 55, 60, field15);  atoms[i].occupancy = (REAL) atof(field15);
	readCol(line, 61, 66, field16);  atoms[i].tempFactor = (REAL) atof(field16);
	readCol(line, 67, 67, field17);
	readCol(line, 68, 70, field18);  fld2i(field18, 3,&(atoms[i].ftNote));
	readCol(line, 71, 72, field19);
	readCol(line, 73, 76, field20);  fld2s(field20, atoms[i].segID);
	readCol(line, 77, 78, field21);  fld2s(field21, atoms[i].element);
	readCol(line, 79, 80, field22);  fld2s(field22, atoms[i].charge);

	/* construct a proper atom name from the info in the ATOM record */
	/* this is not yet perfect */
	strcpy(atoms[i].atomName,atoms[i].atomType);
	strcat(atoms[i].atomName,atoms[i].atomLoc);

/* 	strcpy(atoms[i].iCode,"\0"); */
	i++;
	assert (i < MAXPDB);
    } 

}  /* end of "while (TRUE)" loop */

fclose(file_ref);
nAtom = i;

#ifdef DEBUG
  fprintf(stderr, "REMARK Number of points = %d", nAtom);
  if (pdb_file) fprintf(stderr," in file %s\n",pdb_file);
  else fprintf(stderr," in stdinput\n");
#endif

/* Copy the PDB structure in another array with just enough memory */
pureAtoms = (PDB *) malloc (nAtom*sizeof(PDB));
for (i=0; i<nAtom; i++)  pureAtoms[i] = atoms[i];
*numAtoms = nAtom;
*Atoms = pureAtoms;
free (atoms);

return;
}  /* end of main */


/*************************************************************************************
 * void readCol() - extract a string from columns i1 to i2 of a larger string
 *************************************************************************************/

void readCol(char * line, int i1, int i2, char * field00 )
{
  int i;
  for (i=i1;i<=i2;i++)  field00[i-i1] = line[i-1];
  field00[i-i1] = '\0';
}


/*************************************************************************************
 * void fld2s() - read a string from a string, terminate with end-of-string marker
 *************************************************************************************/

void fld2s(char * field, char * str)
{
  int i;
  i = sscanf(field, "%s", str);
  if  (i < 1) { str[0] = '\0'; }
  /* if  (i < 1) { str[0] = ' '; str[1]='\0'; } */
}

/*************************************************************************************
 * void fld2i() - read integers from a string
 *************************************************************************************/

void fld2i(char * field, int n, int * num)
{
  int i;
  for (i=0; i<n; i++)  {
    if (field[i] != ' ')  { sscanf(field, "%d", num); return; }
  }
  if (i == n)  { *num = 0;  return; }
}


/*************************************************************************************
 *   initVDW.c
 *     read energy parameters 
 *     input : filename
 *     output: atoms.LJ_a, atoms.LJ_b, atoms.LJ_c
 *************************************************************************************/

void initVDW (char* filename, int numAtoms, PDB *atoms)

{
static int New=TRUE;
static int Num;
FILE *ffl;
PDB * anatom;
	
char line[MAXCHARLINE+1];

    int i;
    AtomParam *dict_record;
    AtomParam *pt;

    if (!New) return;


    ffl = fopen(filename, "r");
    if (ffl==NULL ) {
	fprintf(stderr, "Can't open file! %s\n", filename); 
	exit(1);
    }

    New=FALSE;

    anatom = atoms;
    while (TRUE) {
	if (fgets (line, MAXCHARLINE, ffl) == NULL)  break;
	if (EQU(line,"ATOM")) {
	    sscanf (line+66, "%f %f %f %f %s %d",
		&anatom->LJ_c, &anatom->LJ_a, &anatom->LJ_b, &anatom->ms_rad,
		anatom->type, &anatom->back);
	    anatom++;
	}
    }
    fclose (ffl);

} /* end initVDW() */

/* NEW LISTING
ATOM     12  NH2 ARG      1     25.367  12.797 -13.838  1.00 25.73  -0.260   47.56   861.0
ATOM     13 HH22 ARG      0     24.865  12.443 -13.049  0.00  0.00   0.240    0.00     0.0
ATOM     14 HH21 ARG      0     24.897  12.933 -14.710  0.00  0.00   0.240    0.00     0.0
01234567890123456789012345678901234567890123456789012345678901234567890123456789
0         1         2         3         4         5         6         7 
*/

/*
REMARK Number of residues = 58
REMARK number of atoms in output molecule=571
ATOM      1  N   ARG      1     31.758  13.358 -13.673  -0.280   47.56   861.0
ATOM      2  H   ARG      1     31.929  13.460 -14.653   0.280    0.00     0.0
012345678901234567890123456789012345678901234567890123*2345678*2345678*2345678
*/
