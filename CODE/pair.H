/* * * * * * * * * * * * * * * * * * * * * * * * * 
 *    S*I*G*M*A - 1996 program
 * revised 1996        J. Hermans, Univ. N. Carolina
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * *
 *  pair.H: data related to pairlist(s)
 * * * * * * * * * * * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * *
 *  SetLongRange()
 *  SetCutOff2()
 *  SetCutOff()
 *  GetCutOff2()
 *  GetCutOff()
 *  SetNumBoxCutOff()
 *  ShowCutOff()
 *  InitPair()
 * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef pair_h
#define pair_h
int Number_of_planes,Number_of_lines,Num_cubes[3];

/* for each atom */
typedef struct { SHORT link,y_index,z_index; } A_Hash;
A_Hash *AtomHash;

#define Y_INDEX(i) AtomHash[i].y_index /* y-value of the atom */
#define Z_INDEX(i) AtomHash[i].z_index /* z-value of the atom */
#define NEXT(i)   AtomHash[i].link    /* following atom along the line  */

/* for each line of atoms */
typedef struct { SHORT firstatom,lowz,highz; } L_Hash;
L_Hash *LineHash;

#define FIRSTINLINE(ix,iy) LineHash[ix*NBOX_Y+iy].firstatom /* first atom in the line */
#define LOWZ_INLINE(ix,iy) LineHash[ix*NBOX_Y+iy].lowz /* lowest z in the line */
#define HIZ_INLINE(ix,iy)  LineHash[ix*NBOX_Y+iy].highz /* highest z in the line */

#define NBOX_Y Num_cubes[1]

int Boxsize; /* how big is the box in naybor (usually 2 Angstrom) */
#endif
