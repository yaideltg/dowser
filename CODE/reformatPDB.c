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
 *   reformatPDB.c
 *      Read in the original pdb file and a file with residue descriptions
 *      Eliminate nonpolar hydrogens, insert polar hydrogrens
 *      Choose proper chain termini (MUST DO)
 *      Detect residues that are not in the dictionary
 *      Read in additional data to help fix the problems (MUST DO)
 *
 *      Read in a file with geometry and connectivity, and compute missing xyz
 *
 *      Output a new pdb-formatted file
 *
 *   Input:
 *      -pdbin FILENAME - original pdb file 
 *      -atomtypes FILENAME or
 *          $DOWSER/DATA/atomdict.db - dictionary file with residue descriptions
 *      -atomparms FILENAME or
 *          $DOWSER/DATA/atomparms.db - dictionary file with LJ parameters 
 *   Output:
 *      -pdbout FILENAME - new PDB file
 *   Additional optional arguments
 *      -hetero = use ATOM and HETATM records in the pdb file
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dowser.h"

#define SSBONDSQ  5.0 /* if S-S distance^2 is below this, it's a bond */
#define NCBONDSQ  4.0 /* if N-C distance exceeds this, it's a chain break */
#define MAXTYPES 500
#define MAXMERGE 100

extern void  readPDB10(char *, int *, PDB **);
extern void  readDict();
extern int   BackOneBond();
extern char  *getenv();
extern int   ForwardOneBond();
extern void  FindSSBonds ();
extern int   CysSG ();
extern int   Add1Atom();
extern void  FindChainBreaks();
extern int   ScanDictVsPDB();
extern void  MergeTerminus();
extern void  readParam ();
extern void  WriteXYZ ();

FILE *outfile;

/* begin0,end0 = begin and end  in the set that remains when heteroatoms are removed 
   begin,end   = begin and end  in the set containing also the heteroatoms
   newbegin,newend = begin and end in ATOMTYPE */
typedef struct {
  char name[5],terminus[5];
  int begin0,end0,begin,end,newbegin,newend,type,numat,numextra;
  int chainend,index; /* 1 for N-term, 2 for C-term */
} RESIDUE;

ATOMTYPE *atomtypes;

void main(int argc, char *argv[])
{
  PDB	*newatoms, *atoms, *anatom;
  RESIDUE *residues, *aresidue;
  ATOMTYPE *anatomtype, *bnatomtype, *useatomtypes;
  RESTYPE *restypes, *arestype;
  ATOMTYPE mergeatomtype[MAXMERGE];
  RESTYPE mergerestype;
  int *type_key;
  int nType, nRestype;

  int i, j, k, ires, iat;
  int nAtom, nRes;
  int newnAtoms;
  int oldresSeq; char oldchainID; char oldinsert;
  char *oldtype;
  char *what;
  int found;
  char junk [10];
  int pdbin = 0, pdbout =0;
  int error=0;

  int noheteros = TRUE;
  char filetypes[256],fileparms[256];

  /* look for arguments that affect the way dowser is to be used */
  *filetypes = *fileparms = NULLCHAR;
  for (i=1; i<argc; i++) {
      if (EQUAL(argv[i],"-hetero")) noheteros = FALSE;
      if (EQUAL(argv[i],"-atomtypes")) strcpy(filetypes,argv[i+1]);
      if (EQUAL(argv[i],"-atomparms")) strcpy(fileparms,argv[i+1]);
      if (EQUAL(argv[i],"-pdbin")) pdbin=i+1;
      if (EQUAL(argv[i],"-pdbout")) pdbout=i+1;
  }

  if (pdbin*pdbout == 0) {
      fprintf(stderr,"ERROR: reformatPDB needs -pdbin FILENAME -pdbout FILENAME\n");
      exit (1);
  }

  if (noheteros) fprintf(stderr,"* DOWSER - default, no HETATM records will be used\n");
  else {
      fprintf(stderr,"* DOWSER - both ATOM and HETATM records will be used\n%s\n",
      		     "*          except atom 'O' in residue 'HOH'");
  }

  /* process the original PDB file */
  /* the readPDB routine allocates the needed storage for atoms[] */
  readPDB10(argv[pdbin], &nAtom, &atoms);

  if (!(outfile = fopen (argv[pdbout],"w"))) {
    fprintf (stderr,"REMARK ERROR: cannot open output pdb file %s\n",argv[pdbout]);
    exit (1);
  }


/*
ATOM     60  C   ASN    25      33.889  33.642  29.803  1.00 12.34           C  
*/
  /* count the number of residues in the ATOM records of the PB file, based on changes in
	 chain number and/or sequence number */
  anatom = atoms;
  nRes = 0; oldresSeq=0; oldchainID='\0'; oldinsert=' ';
  for (iat=0; iat<nAtom; iat++) {

    if ( !EQUAL(anatom->recdName,"ATOM") &&
         !EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }
    if ( noheteros ) {
	if (EQUAL(anatom->recdName,"HETATM")) {
	    if (!EQUAL(anatom->resName,"HOH")) {
		fprintf(stderr,"* HETERO atom %s %s %s %d%s will be skipped\n",
		anatom->atomType, anatom->resName, anatom->chainID, anatom->resSeq, anatom->iCode);
	    }
	    anatom++; continue;
	}
    }
    else if (EQUAL(anatom->recdName,"HETATM") && EQUAL(anatom->resName,"HOH"))
	{ anatom++; continue; }

    if (oldresSeq != anatom->resSeq || oldchainID != anatom->chainID[0] || oldinsert != anatom->iCode[0])
	{ nRes++ ; oldresSeq = anatom->resSeq; oldchainID = anatom->chainID[0];
	  oldinsert = anatom->iCode[0]; }
    anatom++;
  }
  fprintf (outfile,"REMARK Number of residues = %d\n",nRes);

  /* allocate the structure describing the residues */
  residues = (RESIDUE *) malloc (nRes * sizeof(RESIDUE));

  /* fill the Residue array from the info in PDB */
  anatom = atoms; aresidue = residues;
  nRes = 0; oldresSeq=0; oldchainID='\0';
  for (iat=0; iat<nAtom; iat++) {

    if ( !EQUAL(anatom->recdName,"ATOM") &&
         !EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }
    if ( noheteros ) {
	if (EQUAL(anatom->recdName,"HETATM")) { anatom++; continue; }}
    else if (EQUAL(anatom->recdName,"HETATM") && EQUAL(anatom->resName,"HOH"))
	{ anatom++; continue; }
	  

    /* new residue if residue number or chain-id changes */
    if (oldresSeq != anatom->resSeq || oldchainID != anatom->chainID[0] || oldinsert != anatom->iCode[0])
	{
	    nRes++ ; oldresSeq = anatom->resSeq; oldchainID = anatom->chainID[0];
	    oldinsert = anatom->iCode[0];
	    strcpy(aresidue->name, anatom->resName);
	    aresidue->begin0 = aresidue->end0 = iat; /* init. end for 1-atom residues */
	    aresidue->chainend = aresidue->numextra = aresidue->numat = 0;
	    aresidue->index = nRes;
	    aresidue++;
	}
    else {
	(aresidue-1)->end0=iat;
    }
    anatom++;
  }
  fprintf(stdout,"* Input pdb file contains %d atoms\n",(int)(anatom-atoms));

#ifdef TEST
  fprintf (stdout,"Residues\n");
  for (i=0;i<nRes;i++) fprintf (stdout,"RES: %s %d %d %d %d\n",residues[i].name,
	  residues[i].begin, residues[i].end,
	  residues[i].begin0, residues[i].end0);
#endif

  /* read in the file describing the atoms with charges, connectivity, etc., by residue */
  (void) readDict (filetypes, &nType, &nRestype, &atomtypes, &restypes);

  /* scan the residue types in the molecule vs. those in the dictionary */
  if (ScanDictVsPDB(nRes,residues,nRestype,restypes)) {
      exit (1);
  }

  /* Locate SS bridges and change the residue names to CSS */
  (void) FindSSBonds (atoms,residues,nRes,restypes,nRestype);

  /* Locate chain breaks and mark these in the residues[] by setting terminus */
  (void) FindChainBreaks(atoms,residues,nRes,restypes,atomtypes);

  /* Figure the extra atoms needed for the chain termini */
  for (ires=0;ires<nRes;ires++) {
      if (residues[ires].chainend && *(residues[ires].terminus)) 
	  (void) MergeTerminus
	      (TRUE,atoms,residues+ires,nRestype,restypes,atomtypes,
	      mergeatomtype,&mergerestype);
  }

  /* compute the number of atoms in the new molecule and allocate storage */
  newnAtoms = 0;
  for (ires=0; ires<nRes; ires++) {
	  aresidue = residues+ires;
	  aresidue->begin = newnAtoms;
	  arestype = restypes + aresidue->type;
	  newnAtoms += arestype->end - arestype->begin + 1 + aresidue->numextra;
	  aresidue->newend = newnAtoms-1;
  }
  fprintf (outfile,"REMARK number of atoms in output molecule = %d\n",newnAtoms);
  newatoms = (PDB *) malloc (newnAtoms * sizeof (PDB));

  /* fill the new array of atoms */
  /* look for atom by same name in the input pdb file */
  /* and use its information; else leave some items blank for later
	computation */
  newnAtoms=0;
  for (ires=0; ires<nRes; ires++) {
      aresidue = residues + ires;
      arestype = restypes + aresidue->type;

      /* is there a terminus to add? */
      if (aresidue->chainend && *(aresidue->terminus)) {
	  (void) MergeTerminus
	      (FALSE,atoms,residues+ires,nRestype,restypes,atomtypes,
	      mergeatomtype,&mergerestype);
	  arestype=&mergerestype;
	  useatomtypes=mergeatomtype;

#ifdef TESTM
    anatomtype=mergeatomtype;
    for (i=0;i<mergerestype.numat; i++) {
	j=anatomtype->index;
	bnatomtype=atomtypes+j;
	fprintf(stdout,"%4d %-4s %-4s %-4s %-4s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	  j,bnatomtype->name,bnatomtype->back,bnatomtype->forward,bnatomtype->type,
	  bnatomtype->charge,bnatomtype->bond,bnatomtype->angle,bnatomtype->dihedral);
	anatomtype ++;
    }
#endif
      }
      else useatomtypes=atomtypes;

      for (j=arestype->begin; j<=arestype->end; j++) {
	  anatomtype = useatomtypes + j; 
	  what = anatomtype->name; 
	  found = 0;
	  for (k=aresidue->begin0; k<= aresidue->end0; k++) {
	      strcpy(junk,atoms[k].atomName); 
	      if (EQUAL(what,junk))  { 
		  found = 1; 
		  break; 
	      }
	  }
	  if (found) {
	      newatoms[newnAtoms] = atoms[k];
	      newatoms[newnAtoms].ftNote = TRUE;
	      /* fprintf(stdout,"FOUND %s in %d %s\n",what,i,aresidue->name); */
	  }
	  else {
	      strcpy(newatoms[newnAtoms].atomName,what);
	      strcpy(newatoms[newnAtoms].resName,arestype->resname);
	      newatoms[newnAtoms].ftNote = FALSE;
	      if (what[0]!='h' && what[0]!='H')
		  fprintf(outfile,"REMARK WARNING non-hydrogen NOT FOUND %s in %d %s\n",
		      what,ires+1,aresidue->name);
	  }
	  newatoms[newnAtoms].newResNo = ires;
	  newatoms[newnAtoms].LJ_c     = anatomtype->charge;
	  newatoms[newnAtoms].key_dict = anatomtype->index;
	  strcpy (newatoms[newnAtoms].atomType, anatomtype->type);
	  newnAtoms++;
      }

  } /* end fill new array of atoms */

#ifdef TESTX
  for (i=0;i<newnAtoms;i++) {
	  fprintf(stdout,"NEWATOM %-4s %-4s ",newatoms[i].resName,newatoms[i].atomName);
	  if (newatoms[i].ftNote) fprintf(stdout," OLD ");
	  else fprintf(stdout," NEW ");
	  j=newatoms[i].key_dict;
	  fprintf(stdout," %4d %-4s %-4s %-4s\n",
	      j,atomtypes[j].back,atomtypes[j].forward,atomtypes[j].type);
  }
#endif

    fprintf (stderr, "Add atoms, pass number 1\n");
    if (i = AddAllAtoms (newnAtoms,newatoms,atomtypes, FALSE)) {
	while (TRUE) {
	    fprintf (stderr, "Add atoms, another pass\n");
 	    j = AddAllAtoms (newnAtoms,newatoms,atomtypes, FALSE);
 	    k = AddAllAtoms (newnAtoms,newatoms,atomtypes, TRUE);
	    if (k == 0 || k == i) break;
	    i = k;
	}
	if ( k>0 ) {
	    fprintf (stderr, "REFORMAT: **** ERROR **** Not all atoms have coordinates\n");
	    exit (1);
	}
    }

    /* read the LJ parameters and assign them */
    (void) readParam (fileparms, newnAtoms, newatoms);

    /* Print atoms and coordinates as roughly PDB */
    (void) WriteXYZ (newnAtoms, newatoms, residues);

    exit (error);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function ForwardOneBond()
 *    Find the Forwardchain from the info in "atoms[]" and in "atom types[]"
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int ForwardOneBond(from,atoms,atomtypes,nAtoms)
int from,nAtoms;
PDB *atoms;
ATOMTYPE *atomtypes;

{
int j;
int found=FALSE;
char *forwardname;
ATOMTYPE *anatomtype;

    anatomtype = atomtypes + atoms[from].key_dict;
    forwardname = anatomtype->forward;
	if (EQUAL(forwardname,"NOT")) return (-1);
    for (j=from+1; j<nAtoms; j++) {
	if (EQUAL(forwardname,atoms[j].atomName)) return (j);
    }
    fprintf(stdout,"REMARK ERROR - FORWARDCHAIN: none found with name %s starting from atom %d\n",
	forwardname,from);
    return (-1);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function BackOneBond()
 *    Find the backchain from the info in "atoms[]" and in "atom types[]"
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int BackOneBond(from,atoms,atomtypes)
int from;
PDB *atoms;
ATOMTYPE *atomtypes;

{
int j;
int found=FALSE;
char *backname;
ATOMTYPE *anatomtype;

    anatomtype = atomtypes + atoms[from].key_dict;
    backname = anatomtype->back;
    if (EQUAL(backname,"NOT")) return (-1);
    for (j=from-1; j>=0; j--) {
	if (EQUAL(backname,atoms[j].atomName)) return (j);
    }
    fprintf(stdout,"REMARK ERROR - BACKCHAIN: none found with name %s starting from atom %d\n",
	backname,from);
    return (-1);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  function readDict
 *    input a file with atom types = a dictionary of atoms in residues
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void readDict (filename, numat, numrestyp, atomtypes, restypes)
int * numat, *numrestyp;
ATOMTYPE **atomtypes;
RESTYPE **restypes;
char *filename;

{
FILE *in_file;
ATOMTYPE *anatomtype;
RESTYPE *arestype;
char recdName[10],resname[5];
char term[5],nterm[5],cterm[5];
int iat,i,j;
char line[256];
char *what;

if (*filename == NULLCHAR) {
    what = getenv ("DOWSER");
    if (what) strcpy (filename,what);
    else {
	fprintf (stderr,"REMARK ERROR: must first set environment variable 'DOWSER'\n");
	exit (1);
    }
    strcat (filename,"/DATA/atomdict.db");
}
if (!(in_file = fopen (filename,"r"))) {
    fprintf (stderr,"REMARK ERROR: cannot open file %s\n",filename);
    exit (1);
}

*atomtypes = (ATOMTYPE *) malloc (MAXTYPES * sizeof(ATOMTYPE));
*restypes = (RESTYPE *) malloc (MAXTYPES * sizeof(RESTYPE));

/*
RESIDUE ALA TERM NH3 COO
ATOM ALA  N    C    CA    1.320  114.0  180.0  -0.280    N
*/

*numat=0; *numrestyp=0;
anatomtype = *atomtypes;
arestype = *restypes;
iat=0;
while (1) {
    fgets (line, 100, in_file); if (feof(in_file)) break;
    sscanf (line, "%s %s", recdName, resname);
    if (EQUAL(recdName,"residue")) {
#ifdef TEST
	fprintf(stderr,"%s\n",resname);
#endif
	term[0]=NULLCHAR;
	sscanf (line, "%s %s %s %s %s", recdName, resname, term, nterm, cterm);
	strcpy(arestype->resname,resname);
	arestype->begin=iat;
	if (EQUAL(term,"TERM")) {
	    strcpy(arestype->nterminus,nterm);
	    strcpy(arestype->cterminus,cterm);
	}
	else *(arestype->nterminus) = *(arestype->cterminus) = NULLCHAR;
	(*numrestyp)++; arestype++;
    }
    if (EQUAL(recdName,"atom")) {
	sscanf (line, "%s %s %s %s %s %f %f %f %f %s",
	    recdName,
	    anatomtype->resname, anatomtype->name, anatomtype->back, anatomtype->forward,
	    &anatomtype->bond, &anatomtype->angle, &anatomtype->dihedral,
	    &anatomtype->charge, anatomtype->type);
	anatomtype->index = iat;
	(arestype-1)->end=iat;
	(*numat)++; anatomtype++; iat++;
    }
}
arestype->end=iat;

#ifdef TEST
arestype = *restypes;
for (i=0; i<*numrestyp; i++) {
    fprintf(stdout,"RES: %s %d %d\n",arestype->resname,arestype->begin,arestype->end);
    for (j=arestype->begin; j<=arestype->end; j++) {
	anatomtype = *atomtypes+j;
	fprintf(stdout,"ATOM: %s %s\n", anatomtype->name,anatomtype->type);
    }
    arestype++;
}
#endif
} /* end readDict() */

/*
ATOM    349  OD1 ASN    43      29.384  10.062  -2.624  1.00  4.58      1BPI 485
******+++++ ****?*** +****?   ********++++++++********OOOOOOBBBBBB ???  ****
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  WriteXYZ() : write the new pdb file with LJ parameters and charges
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void WriteXYZ(num,atoms,residues)

int num;
PDB * atoms;
RESIDUE *residues;

{
int i, j, rnum;
PDB *anatom;
char aname[5],achar;
char *what;
int backchain;

    for (i=0;i<num;i++) {
	anatom = atoms +i;
	if (!anatom->ftNote) continue;
	rnum = anatom->newResNo;
	strcpy (aname,anatom->atomName);
	if (strlen(aname) == 4) {
	    achar = aname[3];
	    aname[3] = NULLCHAR;
	}
	else achar = ' ';
	what=atomtypes[anatom->key_dict].back;
	if (EQUAL(what,"NOT")) backchain=0;
	else {
	    for (j=i-1;j>=0;j--) if (EQUAL(atoms[j].atomName,what)) break;
	    if (j==-1) {
		backchain=0;
		fprintf(stderr,"Bad backchain atom %d %s, %s\n",i,anatom->atomName,what);
	    }
	    else backchain = i-j;
	}
	fprintf(outfile, 
	    "ATOM  %5d %c%-4s%-3s %5d    %8.3f%8.3f%8.3f  %4.2f%6.2f %7.3f %7.2f %7.1f %6.2f %-4s %3d\n",
	    i+1, achar, aname, residues[rnum].name,
	    rnum + 1, anatom->XX, anatom->YY, anatom->ZZ,
	    anatom->occupancy, anatom->tempFactor,
	    anatom->LJ_c, anatom->LJ_a, anatom->LJ_b, anatom->ms_rad,
	    atomtypes[anatom->key_dict].type,backchain);
    }
    fprintf (stdout,"* Reformatted and extended pdb file contains %d atoms\n", num);
}
/*
ATOM      2  CA  ACE     1      -2.971  -0.162   1.364
123456789012345678901234567890123456789012345678901234567890
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  readParam() : 
 *     read the LJ parameters and assign them to atoms
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void readParam(filename,nAtoms,atoms)
char *filename;
int nAtoms;
PDB *atoms;

{
PDB *anatom;
char line[256];
char *what;
FILE *ffl;
int eof = FALSE;
char type[10];
REAL a, b, ms_rad;
int i;

    if (*filename == NULLCHAR) {
	what = getenv ("DOWSER");
	if (what) strcpy (filename,what);
	else {
	    fprintf (stderr,"REMARK ERROR: must first set environment variable 'DOWSER'\n");
	    exit (1);
	}
	strcat (filename,"/DATA/atomparms.db");
    }

    if (!(ffl = fopen (filename,"r"))) {
	fprintf (stderr,"REMARK ERROR: cannot open file %s\n",filename);
	exit (1);
    }

    while (TRUE) {
	fgets(line, 100, ffl);
	if (feof(ffl)) break;
	if (!strncmp(line,"TYPE",4)) {
	    sscanf (line+7, "%s", type);
	    sscanf (line+10, "%f", &a);
	    sscanf (line+19, "%f", &b);
	    sscanf (line+29, "%f", &ms_rad);
	    anatom = atoms;
	    for (i=0; i<nAtoms; i++) {
		if (EQU(anatom->atomType,type)) {
		    anatom->LJ_a = a;
		    anatom->LJ_b = b;
		    anatom->ms_rad = ms_rad;
		    anatom->ftNote = 2;
		}
		anatom ++;
	    }
	}
    }
    fclose (ffl);
}
/*
REMARK atomtype LJ-a LJ-b    MS-radius
1234567890123456789012345678901234
TYPE   N      49.36   1300.0  2.40
TYPE   N      49.36   1300.0
012345678901234567890123456789
*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   FindSSBonds ( )
 *     Use a distance criterion 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FindSSBonds (atoms,residues,nRes,restypes,nRestype)
PDB *atoms;
RESIDUE *residues;
RESTYPE *restypes;
int nRes,nRestype;

{
int i,i1,i2,j,k;
REAL xs1[3],xs2[3];
REAL ds, ds2;
static int css_type=-1;
char *what;

if (css_type < 0) {
    for (css_type=0; css_type < nRestype; css_type++ )
	if (EQUAL("CSS",restypes[css_type].resname)) break;
    if (css_type==nRestype) {
	fprintf(stderr,"REMARK ERROR: Residue type CSS not found in the type dictionary\n");
	fprintf(stderr,"REMARK ERROR: Will not be able to handle SS-bonds\n");
    }
}
if (css_type==nRestype) return;

  /* locate disulfide bridges and change residue names */
  for (i1=0;i1<nRes;i1++) {
      if (j=CysSG (atoms,i1,residues,xs1)) {
	  /* look for another CYS residue */
	  for (i2=i1+1; i2<nRes; i2++) {
	      if (k=CysSG (atoms,i2,residues,xs2)) {
		  ds2=0.;
		  for (i=0;i<3;i++) { ds = xs2[i] - xs1[i]; ds2 += ds*ds; }
		  if ( ds2 <= SSBONDSQ ) {
		      /* strcpy(residues[i1].name,"CSS"); */
		      residues[i1].type=css_type;
		      /* strcpy(residues[i2].name,"CSS"); */
		      residues[i2].type=css_type;
fprintf(outfile,"REMARK SS-BOND residues %d and %d, ds2 = %f\n",i1+1,i2+1,ds2);
		      break;
		  }
	      }
	  }
      }
  }
} /* end FindSSBonds() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  int CysSG (atoms,resnum,atoms,xs1)
 *  return atom index of SG if residue is a CYS, else 0
 *  and the coordnates of this atom
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int CysSG (atoms,resnum,residues,xs1)
PDB *atoms;
int resnum; RESIDUE *residues;
REAL *xs1;

{
int j=0;
  if (EQUAL(residues[resnum].name,"CYS")) {
      if (j=LocateInResidue(atoms,resnum,residues,"SG")) {
	  xs1[0] = atoms[j].XX;
	  xs1[1] = atoms[j].YY;
	  xs1[2] = atoms[j].ZZ;
      }
  }
  return (j);
} /* end CysSG() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  int LocateInResidue(atoms,resnum,residues,name)
 *     return index of atom with given name
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int LocateInResidue(atoms,resnum,residues,name)
PDB *atoms;
int resnum; RESIDUE *residues;
char *name;

{
int j;
    for (j=residues[resnum].begin0;j<=residues[resnum].end0;j++) {
      if (EQUAL(atoms[j].atomName,name)) return (j);
    }
    return (0);
} /* end LocateInResidue() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  FindChainBreaks(atoms,residues,nRes,restypes)   
 *     Use a distance criterion 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void FindChainBreaks(atoms,residues,nRes,restypes,atomtypes)
PDB *atoms;
RESIDUE *residues;
int nRes;
RESTYPE *restypes;
ATOMTYPE *atomtypes;

{
int ires,i;
int j1,j2;
REAL dlink, dlink2;
int num_breaks=0;
RESTYPE *arestype;
char *aname,*bname;
int type,begin;

/* chain termination for first and last molecule */
arestype = restypes + residues[0].type;
strcpy(residues[0].terminus,arestype->nterminus);
residues[0].chainend=1;
arestype = restypes + residues[nRes-1].type;
strcpy(residues[nRes-1].terminus,arestype->cterminus);
residues[nRes-1].chainend=2;

for (ires=1;ires<nRes;ires++) {
    /* first atom in the residue according to the dictionary */
    type = residues[ires].type;
    begin = restypes[type].begin;
    aname = atomtypes[begin].name;

    j2 = LocateInResidue(atoms,ires,residues,aname);
    bname = atomtypes[begin].back;
    if (EQUAL(bname,"NOT")) j1=0; /* always a break */
    else j1 = LocateInResidue(atoms,ires-1,residues,"C");

    /* check distance of link */
    dlink2=0;
    if (j1 && j2) {
	dlink = atoms[j1].XX - atoms[j2].XX;
	dlink2 += dlink*dlink;
	dlink = atoms[j1].YY - atoms[j2].YY;
	dlink2 += dlink*dlink;
	dlink = atoms[j1].ZZ - atoms[j2].ZZ;
	dlink2 += dlink*dlink;
	if (dlink2 <= NCBONDSQ) continue; /* it's a bond */
    }

    fprintf (outfile,"REMARK CHAIN BREAK at residue %d, C-N 'bondlength' = %8.2f\n", ires+1, sqrt(dlink2));
    num_breaks++;

    arestype = restypes + residues[ires].type;
    residues[ires].chainend=1;
    strcpy(residues[ires].terminus,arestype->nterminus);
    arestype = restypes + residues[ires-1].type;
    residues[ires-1].chainend=2;
    strcpy(residues[ires-1].terminus,arestype->cterminus);
}

if (!num_breaks) fprintf(outfile,"REMARK NO CHAIN BREAKS\n");

} /* end FindChainBreaks() */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    void MergeTerminus
 *        (countonly,atoms,aresidue,nRestype,restypes,atomtypes,
 *         mergeatomtype,mergerestype)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void MergeTerminus
    (countonly,atoms,aresidue,nRestype,restypes,atomtypes,mergeatomtype,mergerestype)
int countonly;
PDB *atoms;
int nRestype;
RESIDUE *aresidue;
RESTYPE  *restypes, *mergerestype;
ATOMTYPE *atomtypes, *mergeatomtype;

{
char *what;
RESTYPE *arestype, *brestype, *res_restype, *term_restype;
int found, i, j;
int keytype;
ATOMTYPE *anatomtype, *batomtype;

    /* find the terminus residue in the array of residue types */
    what=aresidue->terminus;
    for (i=0;i<nRestype;i++) {
	arestype = restypes+i;
	if (EQUAL(arestype->resname,aresidue->terminus)) term_restype=arestype;
	if (EQUAL(arestype->resname,aresidue->name))     res_restype=arestype;
    }

    if (countonly) {
	for (i=term_restype->begin;i<=term_restype->end;i++) {
	    what=atomtypes[i].name;
	    found=FALSE;
	    for (j=res_restype->begin;j<=res_restype->end;j++) {
		if (EQUAL(atomtypes[j].name,what)) {
		    found = TRUE; break;
		}
	    }
	    if (!found) aresidue->numextra++;
	}
	if (aresidue->numextra)
	fprintf (outfile,"REMARK Need %d extra atoms in residue %s %d\n",
	    aresidue->numextra, aresidue->name, aresidue->index);
    }

    else {
	keytype=0;
	if (aresidue->chainend==1) { /* insert N-terminus */
	    for (i=term_restype->begin;i<=term_restype->end;i++)
		mergeatomtype[keytype++] = atomtypes[i];
	} /* end chain begin */
	arestype = restypes + aresidue->type;
	for (i= arestype->begin; i<=arestype->end; i++) {
	    anatomtype = atomtypes + i;
	    what = anatomtype->name; 
	    for (batomtype=mergeatomtype;batomtype<mergeatomtype+keytype;batomtype++)
		if (EQUAL(batomtype->name,what)) goto SKIP_ME;
	    mergeatomtype[keytype++] = *anatomtype;
	    SKIP_ME: ;
	} /* end residue */
	if (aresidue->chainend==2) {
	    for (i=term_restype->begin;i<=term_restype->end;i++) {
		anatomtype = atomtypes + i;
		what = anatomtype->name; 
		for (batomtype=mergeatomtype;batomtype<mergeatomtype+keytype;batomtype++) {
		    if (EQUAL(batomtype->name,what)) {
			*batomtype = *anatomtype; /* overwrite the specs */
			goto SKIP_ME2;
		    }
		}
		mergeatomtype[keytype++] = *anatomtype;
		SKIP_ME2: ;
	    }
	} /* end chainend */
    }
    strcpy (mergerestype->resname,aresidue->name);
    mergerestype->begin=0;
    mergerestype->end= keytype-1;
    mergerestype->numat= keytype;

    return ;

} /* end MergeTerminus */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    int ScanDictVsPDB(nRes,residues,nRestype,restypes)
 *    compare the sequence in the PDB wioth the available residue types
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int ScanDictVsPDB(nRes,residues,nRestype,restypes)
int nRes,nRestype;
RESIDUE *residues;
RESTYPE *restypes;

{
int i,j;
int error=FALSE;
char *what;

    for (i=0; i<nRes; i++) {
	what = residues[i].name;
	for (j=0; j < nRestype; j++ ) {
	    if (EQUAL(what,restypes[j].resname)) { residues[i].type = j; break; }
	}
	if (j==nRestype) {
	    fprintf(stderr,"REMARK ERROR: Residue type %s not found in the type dictionary\n",what);
	    error=TRUE;
	}
    }
    return (error);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *    int AddAllAtoms()
 *        Compute the coordinates of the missing atoms 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int AddAllAtoms (nAtoms,atoms,atomtypes, final_pass)
int nAtoms;
PDB *atoms;
ATOMTYPE *atomtypes;
int final_pass;

{
int i, fault = 0;
int refatom, backatom, backback;
ATOMTYPE *anatomtype, *bnatomtype, *cnatomtype;

for (i=0;i<nAtoms; i++) {
    if (!atoms[i].ftNote) {
	atoms[i].occupancy = 0.;
	atoms[i].tempFactor = 0.;
	anatomtype = atomtypes + atoms[i].key_dict;
	/* use the single and double backchain for backatom and  backback
		    and the forward chain for refat */
	backatom = BackOneBond (i,atoms,atomtypes);

	if (backatom<0 || final_pass) { /* root atom without coordinates */
	    backatom = ForwardOneBond (i,atoms,atomtypes,nAtoms);
	    if (backatom<0) continue;
	    backback = ForwardOneBond (backatom,atoms,atomtypes,nAtoms);
	    if (backback<0) continue;
	    refatom =  ForwardOneBond (backback,atoms,atomtypes,nAtoms);
	    if (refatom<0) continue;
	    anatomtype =  atomtypes + atoms[refatom].key_dict;
	    bnatomtype =  atomtypes + atoms[backback].key_dict;
	    cnatomtype =  atomtypes + atoms[backatom].key_dict;
	    if (Add1Atom (atoms,backback, backatom, refatom, i,
	      cnatomtype->bond,
	      (PI/180.)*bnatomtype->angle,
	      (PI/180.)* anatomtype->dihedral) ) {
		atoms[i].ftNote = TRUE;
	    }
	    continue;
	}

	backback = BackOneBond (backatom,atoms,atomtypes);
	refatom = ForwardOneBond (backatom,atoms,atomtypes,nAtoms);
	if (refatom > i && atoms[refatom].atomName[0] != 'H') {
	    bnatomtype =  atomtypes + atoms[refatom].key_dict;
	    if (Add1Atom (atoms,backback, backatom, refatom, i,
	      anatomtype->bond, (PI/180.)*anatomtype->angle,
	      (PI/180.)* (anatomtype->dihedral - bnatomtype->dihedral)) )
	      atoms[i].ftNote = TRUE;
	}
	else {
	    /* use 1, 2 and 3 times backchain for backatom, backback and back(backback),
	    last to be used as refat */
	    if (backback>=0) refatom = BackOneBond (backback,atoms,atomtypes);
	    if (backback>=0 && refatom >=0) {
	      if (Add1Atom (atoms,backback, backatom, refatom, i,
		anatomtype->bond, (PI/180.)*anatomtype->angle, (PI/180.)*anatomtype->dihedral) )
		atoms[i].ftNote = TRUE;
	    }
	}
    }
}

for (i=0;i<nAtoms; i++) {
    fault += ( 1 - atoms[i].ftNote );
}
fprintf (stderr, "%d atoms not found\n", fault);
return fault;
} /* end AddAllAtoms() */

char *Test()

{
    char *a;
    a = malloc(100);
    return a;
}
