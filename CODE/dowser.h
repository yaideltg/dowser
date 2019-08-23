#ifndef dowser_h
#define dowser_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#ifndef	DOWSER_H
#define	DOWSER_H
#define WideLine        133 /*  For files that are not exactly 80 colum wide */
#define	MaxRadius	3.2
#define	RaySpacing	0.2
#define	MINIMUM -1000.0
#define	MAXMUM   1000.0
#define MAXSHIFTSQ 16.0 /* limit on displacement in placeWat() */

#define TOL 1.0e-4   /* Tolerance on line searches w/i brent */
#define FTOL 1.0e-4  /* Tolerance on the min. function return in frprmn() */

#define CALL (void)

#define MAXPAIRS 1000         /* Max. # of water-protein pairs */
#define REAL float

#define EQUAL !strcasecmp
#define EQU(a,b)   !strncasecmp(a,b,strlen(b))
#define NULLCHAR '\0'
#define	TRUE	1
#define FALSE	0
#define	PI	3.1415926535

#define CUT_OFF 12         /* cut-off for water NBPairs */

#define A_OW 51.0             /* vdw parameters for water in kcal/mol */
#define B_OW 1623.0
#define C_OW -0.82
#define A_H 0.0
#define B_H 0.0
#define C_H  0.41

#define XX xyz[0]
#define YY xyz[1]
#define ZZ xyz[2]
typedef struct pdbRec
{
    char  recdName[7];    /*       1 -  6 */
    int   serial;         /*       7 - 11 */
    /*char  dummy02[2];     **           12 */
    /*char  name[5];        **      13 - 16 */
    char  atomType[3];
    char  atomLoc[3];
    char  altLoc[2];      /*           17 */
    char  resName[4];     /*      18 - 20 */
    /*char  dummy06[2];     **           21 */
    char  chainID[2];     /*           22 */
    int   resSeq;         /*      23 - 26 */
    char  iCode[2];       /*           27 */
    /*char  dummy10[4];     **      28 - 30 */
    REAL xyz[3];
/*    REAL x;              **      31 - 38 */
/*    REAL y;              **      39 - 46 */
/*    REAL z;              **      47 - 54 */
    REAL occupancy;      /*      55 - 60 */
    REAL tempFactor;     /*      61 - 66 */
    /*char  dummy16[2];     **           67 */
    int   ftNote;         /*      68 - 70 */
    /*char  dummy18[3];     **      71 - 72 */
    char  segID[5];       /*      73 - 76 */
    char  element[3];     /*      77 - 78 */
    char  charge[3];      /*      79 - 80 */

/* the following items are not standard PDB, and are first written by reformatPDB program */
/* subsequent programs may access this information using initVDW() routine */
    char atomName[5];     /* constructed from above */
    int  key_dict;        /* link to additional info */
    int newResNo;         /* number in new order */
    REAL LJ_a, LJ_b;      /* sqrt(LJ coefficients) */
    REAL ms_rad;          /* radius for MS program */
    REAL LJ_c;            /* the atom charge */ 
    char type[5];         /* atom type, cf. atomparms.db */
    int back;             /* backchain */
} PDB;

typedef struct water { PDB Ow; PDB H1; PDB H2; } WATER;

typedef struct {
    char resname[5],name[5],back[5],forward[5],type[5];
    REAL a,b,charge,bond,angle,dihedral,x,y,z;
    int index,backn,forwardn;
} ATOMTYPE;

typedef struct xyz {REAL xyz[3]; } XYZ;

typedef struct ijk{int x; int y; int z;} IJK;

typedef struct sphere { 
    char *atomType; char *resiType; REAL x; REAL y; REAL z; REAL radius;
} Sphere;

/* begin, end = first, last atom in atomtype[] for this residue */
typedef struct {
    char resname[5],nterminus[5],cterminus[5]; int begin,end,numat;
} RESTYPE;


typedef struct
{
    REAL a, b ,c;
} AtomParam;	/*	vdw and charge	*/

#endif
#endif
