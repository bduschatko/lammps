/* ----------------------------------------------------------------------

- Blake Duschatko

- For energy calculations we don't need full system energy
- Energy differences:

    dU = U(N) - U(N+1)
       = sum(intra molecule N) + sum(inter molecular N)
         - sum(intra molecular N) - sum(intra molecular 1)
         - sum(inter molecular N) - sum(inter molecular N+1)
      = - ( sum(intra molecular 1) + sum(inter molecular N+1) )

    Using full energy computes the intra-molecular 1 term

------------------------------------------------------------------------- */

#include "fix_widom.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h" // MPI Communicator
#include "compute.h" // Computes
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h" // Fixes
#include "fix_widom.h"
#include "force.h"
#include "group.h" // Groups
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h" // Math Routines
#include "memory.h"
#include "modify.h"
#include "molecule.h" // Molecule Class
#include "neighbor.h"
#include "pair.h"
#include "random_park.h" // Random Numbers
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

enum { EXCHATOM, EXCHMOL };    // exchmode

/* ----------------------------------------------------------------------
   Initialize Lammps Fix
------------------------------------------------------------------------- */
FixWidomBiased::FixWidomBiased(LAMMPS *lmp, int narg, char **arg) :
    FixWidom::FixWidom(LAMMPS *lmp, int narg, char **arg){

      if (triclinic)
        error->all(FLERR, "Fix widom/biased does not support triclinic cells");  
      if (exchmode != EXCHMOL)
        error->all(FLERR, "Fix widom/biased only supports molecule insertions");

    };

void FixWidomBiased::init()
{}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixWidomBiased::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix widom/biased command");

  // Set up volume biasing and remove from args for parent class
  int iarg, iarg_reduced = 0;
  char reduced_arg[narg-2];
  while (iarg < narg) {
    if (strcmp(arg[iarg],"volume") == 0) {
      if (iarg+1 > narg)
        error->all(FLERR, "Illegal fix widom/biased command");
      read_volumes(arg[iarg+1])
      iarg += 2;
    }
    else
      reduced_arg[iarg_reduced] = arg[iarg]
      iarg_reduced += 1;
  }

  // call parent without volume info 
  FixWidom::options(narg-2, reduced_arg);

}

void FixWidomBiased::read_volumes(char *filename){
  int me = comm->me;
  char line[MAXLINE];
  FILE *fptr;
  char *ptr;

  // Read file on process 0
  if (me = 0){
    fpt = utils::open_potential(filename, lmp, nullptr) // Where is open potential ?
    if (fptr == NULL) {
      char str[128];
      snprintf(str, 128, "Cannot open volumes file %s", filename)
      error->one(FLERR, str);
    }
  }

  int nvols, ncols; 

  if (me == 0){

    // Read up to MAXLINES from file stream, until a new line is reached
    fgets(line, MAXLINE, fptr); 

    // Get number of volumes and column size
    sccanf(line, "%i %i", &nvols, &ncols)

    // Create buffer for file
    int buffer_size = nvols*ncols;
    memory->create(volume_buffer, buffer_size, "fix:buffer")

    for (int i=0; i < nvols; i++){
      fgets(line, MAXLINE, fptr); 
      ptr = strtok(line, " \t\n\r\f");
      volume_buffer[i] = atof(ptr)
      while ((ptr = strtok(NULL, " \t\n\r\f")))
        volume_buffer[i++] = atof(ptr);
    }

  }

  // Broadcast buffer
  MPI_Bcast(&volumes, , MPI_INT, 0, world)

  // Locally reshape buffer 
  volumes = Eigen::MatrixXd::Zero(nvols, ncols);
  for (int i = 0; i < nvols; i++){
      for (int j = 0; j < ncols; j++){
          volumes[i,j] = volume_buffer[i*ncols + j]
      }
  }

}

/* ----------------------------------------------------------------------
   Memory Management 
------------------------------------------------------------------------- */

FixWidomBiased::~FixWidomBiased():
  FixWidom::~FixWidom()
{
  memory->destroy(volume_buffer); 
}

/* ----------------------------------------------------------------------
   Insertion Routines
------------------------------------------------------------------------- */
  
  // Allegro is special case
  // This would normally compute pair interactions one at a time 
  // Allegro better if we parallel

  void FixWidomBiased::attempt_molecule_insertion(){
    /*
     * Perform a molecule insertion and compute pair energy
     * Single point energy computations, no need to create atoms
     */

    // Randomly Sample Volume and Get Center
    int volid; // ID
    double size = volumes[volid][4];

    double com_coord[3];
    com_coord[0] = volumes[volid][0] + random_equal->uniform() * size;
    com_coord[1] = volumes[volid][1] + random_equal->uniform() * size;
    com_coord[2] = volumes[volid][2] + random_equal->uniform() * size;

    // Create a Rotation Matrix for Molecule Orientation
    // --------------------------------------------
    double r[3],rotmat[3][3],quat[4];
    double rsq = 1.1;
    while (rsq > 1.0) {
      r[0] = 2.0*random_equal->uniform() - 1.0;
      r[1] = 2.0*random_equal->uniform() - 1.0;
      r[2] = 2.0*random_equal->uniform() - 1.0;
      rsq = MathExtra::dot3(r, r);
    }
    double theta = random_equal->uniform() * MY_2PI;
    MathExtra::norm3(r);
    MathExtra::axisangle_to_quat(r,theta,quat);
    MathExtra::quat_to_mat(quat,rotmat);
    // --------------------------------------------

    // Place Each Atom in Molecule
    for (int i = 0; i < natoms_per_molecule; i++) {

      // Rotate Molecule Orientation Vector
      MathExtra::matvec(rotmat,onemol->x[i],xtmp);

      double xtmp[3];
      xtmp[0] += com_coord[0];
      xtmp[1] += com_coord[1];
      xtmp[2] += com_coord[2];

      // Remap into Domain
      imageint imagetmp = imagezero;
      domain->remap(xtmp,imagetmp);
      if (!domain->inside(xtmp))
        error->one(FLERR,"Fix widom put atom outside box");

      // Only Create Atom in the Right Subdomain
      int proc_flag = 0;

      if (proc_flag) {
        int ii = -1;
        if (charge_flag) {
          ii = atom->nlocal + atom->nghost;
          if (ii >= atom->nmax) atom->avec->grow(0);
          atom->q[ii] = charge;
        }
        double insertion_energy = energy(ii,nwidom_type,-1,coord);
        double inst_chem_pot = exp(-insertion_energy*beta);
        double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
        ave_widom_chemical_potential += incr_chem_pot / (imove + 1);
      }

    }
  }


  void FixWidomBiased::attempt_molecule_insertion_full(){
    /*
     * Perform a molecule insertion and compute full energy
     */

    double lamda[3];

    // Create new molecule tag
    tagint maxmol = 0;

    // Determine current max molecule tag in system
    // Loop over atoms local to process 
    for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);

    // Collect maxmol tags from all processes 
    tagint maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

    // Increment maxmol tag 
    maxmol_all++;
    if (maxmol_all >= MAXTAGINT)
      error->all(FLERR,"Fix widom ran out of available molecule IDs");
    int insertion_molecule = maxmol_all;

    // Repeat checks for atom number tags
    tagint maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
    tagint maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

    // Randomly Sample Volume and Get Center
    int volid; // ID
    double size = volumes[volid][4];

    double com_coord[3];
    com_coord[0] = volumes[volid][0] + random_equal->uniform() * size;
    com_coord[1] = volumes[volid][1] + random_equal->uniform() * size;
    com_coord[2] = volumes[volid][2] + random_equal->uniform() * size;

    // Create a Rotation Matrix for Molecule Orientation
    // --------------------------------------------
    double r[3],rotmat[3][3],quat[4];
    double rsq = 1.1;
    while (rsq > 1.0) {
      r[0] = 2.0*random_equal->uniform() - 1.0;
      r[1] = 2.0*random_equal->uniform() - 1.0;
      r[2] = 2.0*random_equal->uniform() - 1.0;
      rsq = MathExtra::dot3(r, r);
    }
    double theta = random_equal->uniform() * MY_2PI;
    MathExtra::norm3(r);
    MathExtra::axisangle_to_quat(r,theta,quat);
    MathExtra::quat_to_mat(quat,rotmat);
    // --------------------------------------------

    // Place Each Atom in Molecule
    for (int i = 0; i < natoms_per_molecule; i++) {

      // Rotate Molecule Orientation Vector
      MathExtra::matvec(rotmat,onemol->x[i],xtmp);

      double xtmp[3];
      xtmp[0] += com_coord[0];
      xtmp[1] += com_coord[1];
      xtmp[2] += com_coord[2];

      // Remap into Domain
      imageint imagetmp = imagezero;
      domain->remap(xtmp,imagetmp);
      if (!domain->inside(xtmp))
        error->one(FLERR,"Fix widom put atom outside box");

      // Only Create Atom in the Right Subdomain
      int proc_flag = 0;
      if (triclinic == 0) {
        if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
            xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
            xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) 
            proc_flag = 1; // Atom DOES belong to subdomain of this process
      } else {
        domain->x2lamda(xtmp,lamda);
        if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
            lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
            lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
      }

        // If the atom belongs to this subdomain, add it
        if (proc_flag) {

          atom->avec->create_atom(onemol->type[i],xtmp);
          int m = atom->nlocal - 1;

          atom->image[m] = imagetmp;
          atom->molecule[m] = insertion_molecule;
          if (maxtag_all+i+1 >= MAXTAGINT)
            error->all(FLERR,"Fix widom ran out of available atom IDs");
          atom->tag[m] = maxtag_all + i + 1;
          atom->v[m][0] = 0;
          atom->v[m][1] = 0;
          atom->v[m][2] = 0;

          atom->add_molecule_atom(onemol,i,m,maxtag_all);
          modify->create_attribute(m);
      }

    }

    // Update Global Attributes
    atom->natoms += natoms_per_molecule;
    if (atom->natoms < 0) error->all(FLERR,"Too many total atoms");
    atom->nbonds += onemol->nbonds;
    atom->nangles += onemol->nangles;
    atom->ndihedrals += onemol->ndihedrals;
    atom->nimpropers += onemol->nimpropers;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

    // energy_after corrected by energy_intra
    // WARNING do NOT include intra-molecular energies
    double insertion_energy = (energy_full() -energy_intra) - energy_stored;
    double inst_chem_pot = exp(-insertion_energy*beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    // Remove atoms when done
    atom->nbonds -= onemol->nbonds;
    atom->nangles -= onemol->nangles;
    atom->ndihedrals -= onemol->ndihedrals;
    atom->nimpropers -= onemol->nimpropers;
    atom->natoms -= natoms_per_molecule;

    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == insertion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

  };


/* ----------------------------------------------------------------------
   Energy Routines
------------------------------------------------------------------------- */
  
double FixWidomBiased::energy(int i, int itype, tagint imolecule, double *coord)
{
  if strmtch(lmp->pair_style, "allegro")
    error->all(FLERR, "Allegro not yet supported for Widom insertion");
  else if strmtch(lmp->pair_style, "allegro64")
    error->all(FLERR, "Allegro not yet supported for Widom insertion");
  else
    return FixWidom::energy(int i, int itype, tagint imolecule, double *coord);

}

}
