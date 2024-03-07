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

#include "atom.h"
#include "atom_vec.h"
#include "comm.h" // MPI Communicator
#include "domain.h"
#include "error.h"
#include "fix.h" // Fixes
#include "fix_widom.h"
#include "fix_widom_biased.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h" // Math Routines
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "pair.h"
#include "random_park.h" // Random Numbers
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>
#include <stdio.h>
#include <iostream> // debugging

// check if allegro is built
#ifdef LMP_PAIR_ALLEGRO_H
#include "allegro.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

#define MAXLINE 1024

enum { EXCHATOM, EXCHMOL };    // exchmode

/* ----------------------------------------------------------------------
   Functions to prep inputs for Parent Class
------------------------------------------------------------------------- */

int reduce_array_size(LAMMPS *lmp, int narg, char **arg){

  int iarg = 0;

  while (iarg < narg){
    if (strcmp(arg[iarg],"volume") == 0)
      break;

    else 
      iarg += 1;
  }

  return iarg;

}

char** reduce_array(LAMMPS *lmp, int narg, char **arg){

  int iarg = 0;

  while (iarg < narg){
    if (strcmp(arg[iarg],"volume") == 0)
      return &arg[0];

    else 
      iarg += 1;
  }

  return arg;

}

bool pair_is_allegro(LAMMPS *lmp){

  std::vector<std::string> fixes = {"allegro", 
                                    "allegro3232",
                                    "allegro3264",
                                    "allegro6432", 
                                    "allegro6464"};

  for (int iMod = 0; iMod < lmp->modify->nfix; iMod++){
    if (std::cout(fixes.begin(), fixes.end(), lmp->modify->fix[iMod]->style) > 0)
      return true;
	}
  return false; 
}

// ---------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Initialize Lammps Fix
------------------------------------------------------------------------- */

FixWidomBiased::FixWidomBiased(LAMMPS *lmp, int narg, char **arg):
    FixWidom::FixWidom(lmp, reduce_array_size(lmp, narg, arg), 
                            reduce_array(lmp, narg, arg))
    {

      // Set up volume biasing and remove from args for parent class
      this->options(narg, arg);

      if (triclinic)
        error->all(FLERR, "Fix widom/biased does not support triclinic cells");  
      if (exchmode != EXCHMOL)
        error->all(FLERR, "Fix widom/biased only supports molecule insertions");

      //full_flag = true;
      //error->warning(FLERR, "Using full energy computations, required by Allegro");

      if (narg < 0) error->all(FLERR,"Illegal fix widom/biased command");

    };

FixWidomBiased::~FixWidomBiased(){
    memory->destroy(volume_buffer);
}

void FixWidomBiased::init()
{
    if (pair_is_allegro(this->lmp))
      allegro_flag = true;

    this->FixWidom::init(); // call this child's parent
}

// does NOT override parent class, this gets its own call 
void FixWidomBiased::options(int narg, char **arg)
{ // char** is a pointer to an array of pointers 
  // arg[0] is a pointer 

  if (narg < 0) error->all(FLERR,"Illegal fix widom/biased command");

  // Set up volume biasing and remove from args for parent class
  int iarg = 0;
 
  while (iarg < narg) {
    // how to determine volumes
    if (strcmp(arg[iarg],"volume") == 0) {

      if (iarg+1 > narg)
        error->all(FLERR, "Illegal fix widom/biased command"); // expect info 

      // read from file
      if (strcmp(arg[iarg+1],"file") == 0){
        if (iarg+2 > narg)
          error->all(FLERR, "Illegal fix widom/biased command"); // expect a file name 
        read_volumes(arg[iarg+1]);
        dynamic_volume = false;
        iarg += 3;
      }

      // define resolution in line
      else if ( (strcmp(arg[iarg+1],"spacing") == 0) || (strcmp(arg[iarg+1],"distance") == 0) ){

        if (strcmp(arg[iarg+1],"spacing") == 0){
          if (iarg+2 > narg)
            error->all(FLERR, "Illegal fix widom/biased command"); // expect a value 
          spacing = utils::numeric(FLERR,arg[iarg+2],false,lmp);

          // check if also distance 
          if ((iarg+3 < narg) && strcmp(arg[iarg+3],"distance") == 0){
              if (iarg+4 > narg)
                error->all(FLERR, "Illegal fix widom/biased command"); // expect a value 
              solvent_buffer = utils::numeric(FLERR,arg[iarg+4],false,lmp);
              iarg += 5;
          }
          else
            iarg += 3;

        }
        else { // distance first 
          if (iarg+2 > narg)
            error->all(FLERR, "Illegal fix widom/biased command");
          solvent_buffer = utils::numeric(FLERR,arg[iarg+2],false,lmp);

          // check if also distance 
          if ((iarg+3 < narg) && strcmp(arg[iarg+3],"spacing") == 0){
              if (iarg+4 > narg)
                error->all(FLERR, "Illegal fix widom/biased command"); // expect a value 
              spacing = utils::numeric(FLERR,arg[iarg+4],false,lmp);
              iarg += 5;
          }
          else
            iarg += 3;

        }

        dynamic_volume = true;
      }

      else
        error->all(FLERR, "Illegal fix widom/biased command");

    }

    // rotations resolution 
    else if (strcmp(arg[iarg],"rot") == 0){
      if (iarg+1 > narg)
        error->all(FLERR, "Illegal fix widom/biased command");
      nphi = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }

    // allow max volume evaluations up to NVOLS
    else if (strcmp(arg[iarg],"max") == 0){
      if (iarg+1 > narg)
        error->all(FLERR, "Illegal fix widom/biased command");
      max_inserts = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
    else
      iarg += 1;
  }

}

void FixWidomBiased::read_volumes(char *filename){

  int me = comm->me;
  char line[MAXLINE];
  FILE *fptr;
  char *ptr;

  // Read file on process 0
  if (me == 0){
    fptr = utils::open_potential(filename, lmp, nullptr); // Where is open potential ?
    if (fptr == NULL) {
      char str[128];
      snprintf(str, 128, "Cannot open volumes file %s", filename);
      error->one(FLERR, str);
    }
  }

  if (me == 0){

    // Read up to MAXLINES from file stream, until a new line is reached
    fgets(line, MAXLINE, fptr); 

    // Get number of volumes
    sscanf(line, "%i", &nvols);

    // Create buffer for file
    int buffer_size = nvols;

    // need to allocate memory of the right size for broadcasting the adress
    memory->create(volume_buffer, buffer_size, "fix:volume_buffer");

    for (int i=0; i < nvols; i++){
      fgets(line, MAXLINE, fptr); 
      ptr = strtok(line, " \t\n\r\f");
      volume_buffer[i] = atof(ptr);
      while ((ptr = strtok(NULL, " \t\n\r\f")))
        volume_buffer[i++] = atof(ptr);
    }
  }

  // Broadcast buffer
  MPI_Bcast(&volume_buffer, 3*nvols, MPI_INT, 0, world);
  MPI_Bcast(&nvols, 1, MPI_INT, 0, world);

  // Locally reshape buffer 
  for (int i = 0; i < nvols; i++){          
    //volumes[i].resize(3);
      volumes.push_back({0., 0., 0.});
      for (int j = 0; j < 3; j++){
          volumes[i][j] = volume_buffer[i*3 + j];
      }
  }
}

void FixWidomBiased::pre_exchange(){
  /*
  Compute volumes dynamically 
  */

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  std::cout << "STEPS " << next_reneighbor << " " << update->ntimestep << "\n";

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  if (triclinic) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  if (dynamic_volume){

    // get this subdomain grid
    std::vector<std::vector<double>> grid;

    // number grid points
    int ngrid_x = int((subhi[0] - sublo[0])/spacing);
    int ngrid_y = int((subhi[1] - sublo[1])/spacing);
    int ngrid_z = int((subhi[2] - sublo[2])/spacing);

    int grid_size = ngrid_x * ngrid_y * ngrid_z;

    for (int i = 0; i < ngrid_x; i++){
      for (int j = 0; j < ngrid_y; j++){
        for (int k = 0; k < ngrid_z; k++){
          double x = sublo[0] + i*spacing;
          double y = sublo[1] + j*spacing;
          double z = sublo[2] + k*spacing;
          grid.push_back({x,y,z});
        }
      }
    }

    // compute valid grid points
    // Do NOT need to check ghost atoms on other processors 
    // to avoid repeat grid points
    std::vector<double> myVolumes;
    for (int g = 0; g < grid_size; g++){
      for (int natom = 0; natom < (atom->nlocal); natom++){

          std::vector<double> tmp;
          tmp = {atom->x[natom][0]-grid[g][0], 
                  atom->x[natom][1]-grid[g][1], 
                  atom->x[natom][2]-grid[g][2]};

          // ----------------------------------------------
          // periodic boundaries, this can be more efficient though
          while (tmp[0] < -(xhi-xlo)/2.){
            tmp[0] += (xhi-xlo);
          }
          while (tmp[1] < -(yhi-ylo)/2.){
            tmp[1] += (yhi-xlo);
          }
          while (tmp[2] < -(zhi-zlo)/2.){
            tmp[2] += (zhi-zlo);
          }

          while (tmp[0] > (xhi-xlo)/2.){
            tmp[0] -= (xhi-xlo);
          }

          while (tmp[1] > (yhi-ylo)/2.){
            tmp[1] -= (yhi-ylo);
          }

          while (tmp[2] > (zhi-zlo)/2.){
            tmp[2] -= (zhi-zlo);
          }

          // ----------------------------------------------

          // add if in free volume
          double d = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
          if (d > solvent_buffer)
            myVolumes.push_back(tmp[0]);
            myVolumes.push_back(tmp[1]);
            myVolumes.push_back(tmp[2]);
      }
    }

    // send local dimension to all ranks
    int local_dim = myVolumes.size();
    int global_dim;
    MPI_Allreduce(&local_dim,&global_dim,1,MPI_INT,MPI_SUM,world);

    // Receive all volumes if rank 0
    int me = comm->me;
    int nranks;
    //MPI_Comm_rank(world, &me);
    MPI_Comm_size(world, &nranks);

    if (me == 0){

      int buffer_size = global_dim;
      memory->create(volume_buffer, global_dim, "fix:volume_buffer"); 

      int ibuffer = 0;
      for (int myV = 0; myV < local_dim; myV++){
        volume_buffer[ibuffer] = myVolumes[myV];
        ibuffer += 1;
      }     
      
      for (int r = 1; r < nranks; r++){
        int rank_dim;
        MPI_Status status;
        //auto rank_volumes = new std::vector<double>;
        std::vector<double> rank_volumes;
        MPI_Recv(&rank_dim, 1, MPI_INT, r, 1, world, &status); // receive this ranks dimension
        MPI_Recv(&rank_volumes, rank_dim, MPI_DOUBLE, r, 2, world, &status); // receive this ranks data

        for (int myV = 0; myV < rank_dim; myV++){
          volume_buffer[ibuffer] = rank_volumes[myV];
          ibuffer += 1;
        }    

      }
    }

    // send volumes if not rank 0
    else {
      for (int r = 1; r < nranks; r++){
        MPI_Send(&local_dim, 1, MPI_INT, 0, 1, world);
        MPI_Send(&myVolumes, local_dim, MPI_DOUBLE, 0, 2, world);
      }
    }
    
    // Broadcast buffer
    MPI_Bcast(&volume_buffer, global_dim, MPI_DOUBLE, 0, world);

    // Reformat volumes on all ranks 
    int row = 0;
    while (3*row < global_dim){
      volumes.push_back({volume_buffer[3*row], volume_buffer[3*row+1], volume_buffer[3*row+2]});
      row += 1;
    }

    nvols = volumes.size();

  }

  ave_widom_chemical_potential = 0.0;

  if (region) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  if (full_flag) {
    
    energy_stored = energy_full();

    if (exchmode == EXCHATOM) {
      attempt_atomic_insertion_full();
    } else {
      attempt_molecule_insertion_full();
    }

    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  } else {

    if (exchmode == EXCHATOM) {
      attempt_atomic_insertion();
    } else {
      attempt_molecule_insertion();
    }

  }
  next_reneighbor = update->ntimestep + nevery;

}

/* ----------------------------------------------------------------------
   Insertion Routines
------------------------------------------------------------------------- */
  
  // Allegro is special case
  // This would normally compute pair interactions one at a time 
  // Allegro better if we parallel

void FixWidomBiased::attempt_atomic_insertion(){
    error->all(FLERR, "fix widom/biased does not support atomic insertion");
}

void FixWidomBiased::attempt_molecule_insertion(){
   /*
     * Perform a molecule insertion and compute pair energy
     * Single point energy computations, no need to create atoms
     */

    // Loop through volumes
    // TO DO - OPTIMIZE OPEN MP?
    
    // Make a vector of insertion values and then add them at the end to use openmp, try it
    // a way to do something similar with Allegro?

    int loop_limit = (nvols < max_inserts) ? nvols : max_inserts;

    int moves = 0;
    for (int vID = 0; vID < loop_limit; vID++)
    {

      //double com_coord[3];
      std::vector<double> com_coord;
      com_coord = volumes[vID];

      for (int iphi = 0; iphi < nphi; iphi++)
      {
        // Create a Rotation Matrix for Molecule Orientation

        // axis

        // MAKE FASTER
        double r[3];
        double rsq = 1.1;
        while (rsq > 1.0) {
          r[0] = 2.0*random_equal->uniform() - 1.0;
          r[1] = 2.0*random_equal->uniform() - 1.0;
          r[2] = 2.0*random_equal->uniform() - 1.0;
          rsq = MathExtra::dot3(r, r);
        }

        // relative angle
        for (int itheta = 0; itheta < nphi; itheta++)
        {

          double theta = random_equal->uniform() * MY_2PI;

          double lamda[3];
          double rotmat[3][3],quat[4];
          MathExtra::norm3(r);
          MathExtra::axisangle_to_quat(r,theta,quat);
          MathExtra::quat_to_mat(quat,rotmat);

          // Place Each Atom in Molecule
          auto procflag = new bool[natoms_per_molecule];
          double insertion_energy = 0.0;
          for (int i = 0; i < natoms_per_molecule; i++) {
              
            double xtmp[3];

            // Rotate Molecule Orientation Vector
            MathExtra::matvec(rotmat,onemol->x[i],xtmp);

            xtmp[0] += com_coord[0];
            xtmp[1] += com_coord[1];
            xtmp[2] += com_coord[2];

            // Remap into Universal Domain
            imageint imagetmp = imagezero;
            domain->remap(xtmp,imagetmp);
            if (!domain->inside(xtmp))
              error->one(FLERR,"Fix widom put atom outside box");

            // check if atom belongs to this MPI process/rank (is in sub domain)
            procflag[i] = false;
            if (triclinic == 0) {
              if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
                  xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
                  xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) procflag[i] = true;
            } else {
              domain->x2lamda(xtmp,lamda);
              if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
                  lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
                  lamda[2] >= sublo[2] && lamda[2] < subhi[2]) procflag[i] = true;
            }

            // --------------------------------------------
            // Compute Insertion Energy
            if (procflag[i]){
              int ii = -1;
              if (onemol->qflag == 1) {
                ii = atom->nlocal + atom->nghost;
                if (ii >= atom->nmax) atom->avec->grow(0);
                atom->q[ii] = onemol->q[i];

              }       
            insertion_energy = energy(ii,onemol->type[i],-1,xtmp);
            }
            // --------------------------------------------
          }

          // Collect insertion energies from all processes 
          double insertion_energy_sum = 0.0;
          MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,
                        MPI_DOUBLE,MPI_SUM,world);

          // the insertion_energy_sum is the variable with the energy of inserting one molecule
          double inst_chem_pot = exp(-insertion_energy_sum*beta);
          double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
          ave_widom_chemical_potential += incr_chem_pot / (moves + 1);

          moves += 1; 

          delete[] procflag;

        }
      }
    }
}

void FixWidomBiased::attempt_atomic_insertion_full(){
    error->all(FLERR, "fix widom/biased does not support atomic insertion");
}

void FixWidomBiased::attempt_molecule_insertion_full(){

   /*
     * Perform a molecule insertion and compute TOTAL energy
     */

    // -----------------------------------------------------
    // Set up new molecule 

    tagint maxmol = 0;
    for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
    tagint maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    maxmol_all++;
    if (maxmol_all >= MAXTAGINT)
      error->all(FLERR,"Fix widom ran out of available molecule IDs");
    int insertion_molecule = maxmol_all;

    tagint maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
    tagint maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

    // -----------------------------------------------------


    int loop_limit = (nvols < max_inserts) ? nvols : max_inserts;

    int moves = 0;
    for (int vID = 0; vID < loop_limit; vID++)
    {

      std::vector<double> com_coord;
      com_coord = volumes[vID];

      for (int iphi = 0; iphi < nphi; iphi++)
      {
        // Create a Rotation Matrix for Molecule Orientation

        // axis

        double r[3];
        double rsq = 1.1;
        while (rsq > 1.0) {
          r[0] = 2.0*random_equal->uniform() - 1.0;
          r[1] = 2.0*random_equal->uniform() - 1.0;
          r[2] = 2.0*random_equal->uniform() - 1.0;
          rsq = MathExtra::dot3(r, r);
        }

        // relative angle
        for (int itheta = 0; itheta < nphi; itheta++)
        {

          double theta = random_equal->uniform() * MY_2PI;

          double lamda[3];
          double rotmat[3][3],quat[4];
          MathExtra::norm3(r);
          MathExtra::axisangle_to_quat(r,theta,quat);
          MathExtra::quat_to_mat(quat,rotmat);

          // Place Each Atom in Molecule
          auto procflag = new bool[natoms_per_molecule];
          for (int i = 0; i < natoms_per_molecule; i++) {
              
            double xtmp[3];

            // Rotate Molecule Orientation Vector
            MathExtra::matvec(rotmat,onemol->x[i],xtmp);

            xtmp[0] += com_coord[0];
            xtmp[1] += com_coord[1];
            xtmp[2] += com_coord[2];

            // Remap into Universal Domain
            imageint imagetmp = imagezero;
            domain->remap(xtmp,imagetmp);
            if (!domain->inside(xtmp))
              error->one(FLERR,"Fix widom put atom outside box");

            // check if atom belongs to this MPI process/rank (is in sub domain)
            procflag[i] = false;
            if (triclinic == 0) {
              if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
                  xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
                  xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) procflag[i] = true;
            } else {
              domain->x2lamda(xtmp,lamda);
              if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
                  lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
                  lamda[2] >= sublo[2] && lamda[2] < subhi[2]) procflag[i] = true;
            }

            // --------------------------------------------
            // Add Atom Attributes
            if (procflag[i]){

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
            // --------------------------------------------
          }

          // --------------------------------------------
          // Insertion energy after all atoms added

          // finish molecule set up 
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

          double insertion_energy = (energy_full() -energy_intra) - energy_stored;
          double inst_chem_pot = exp(-insertion_energy*beta);
          double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
          ave_widom_chemical_potential += incr_chem_pot / (moves + 1);

          moves += 1; 

          // --------------------------------------------
          // Clean up, delete molecule 

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

          delete[] procflag;

        }
      }
    }
}

void FixWidomBiased::attempt_molecule_insertion_allegro(){
  #ifndef LMP_PAIR_ALLEGRO_H
    error->all(FLERR,"LAMMPS has not been built with Allegro");
  #else




  #endif
}
