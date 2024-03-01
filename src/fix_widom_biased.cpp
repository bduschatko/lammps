// clang-format off

/* ----------------------------------------------------------------------
   Contributing author: Evangelos Voyiatzis (Royal DSM)
------------------------------------------------------------------------- */

#include "fix_widom.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

enum { EXCHATOM, EXCHMOL };    // exchmode

/* ---------------------------------------------------------------------- */

FixWidom::FixWidomBiased(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), full_flag(false), molcoords(nullptr),
    molq(nullptr), molimage(nullptr), random_equal(nullptr), c_pe(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix widom", error);

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR, "Fix widom does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  triclinic = domain->triclinic;

  // required args

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  ninsertions = utils::inumeric(FLERR,arg[4],false,lmp);
  nwidom_type = utils::inumeric(FLERR,arg[5],false,lmp);
  seed = utils::inumeric(FLERR,arg[6],false,lmp);
  insertion_temperature = utils::numeric(FLERR,arg[7],false,lmp);

  if (nevery <= 0) error->all(FLERR,"Invalid fix widom every argument: {}", nevery);
  if (ninsertions < 0) error->all(FLERR,"Invalid fix widom insertions argument: {}", ninsertions);
  if (seed <= 0) error->all(FLERR,"Invalid fix widom seed argument: {}", seed);
  if (insertion_temperature < 0.0)
    error->all(FLERR,"Invalid fix widom temperature argument: {}", insertion_temperature);

  // read options from end of input line

  options(narg-8,&arg[8]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi = region_zlo = region_zhi = 0.0;
  if (region) {
    if (region->bboxflag == 0)
      error->all(FLERR,"Fix widom region {} does not support a bounding box", region->id);
    if (region->dynamic_check())
      error->all(FLERR,"Fix widom region {} cannot be dynamic", region->id);

    region_xlo = region->extent_xlo;
    region_xhi = region->extent_xhi;
    region_ylo = region->extent_ylo;
    region_yhi = region->extent_yhi;
    region_zlo = region->extent_zlo;
    region_zhi = region->extent_zhi;

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (region->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo) * (region_yhi - region_ylo)
      * (region_zhi - region_zlo);

    region_volume = max_region_volume * static_cast<double>(inside) / static_cast<double>(attempts);
  }

  // error check and further setup for exchmode = EXCHMOL

  if (exchmode == EXCHMOL) {
    if (onemol->xflag == 0)
      error->all(FLERR,"Fix widom molecule {} must have coordinates", onemol->id);
    if (onemol->typeflag == 0)
      error->all(FLERR,"Fix widom molecule {} must have atom types", onemol->id);
    if (nwidom_type != 0)
      error->all(FLERR,"Atom type must be zero in fix widom mol command");
    if (onemol->qflag == 1 && atom->q == nullptr)
      error->all(FLERR,"Fix widom molecule {} has charges, but atom style does not", onemol->id);

    onemol->check_attributes();
  }

  if (charge_flag && atom->q == nullptr)
    error->all(FLERR,"Fix widom atom has charge, but atom style does not");

  // setup of array of coordinates for molecule insertion

  if (exchmode == EXCHATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemol->natoms;
  nmaxmolatoms = natoms_per_molecule;
  grow_molecule_arrays(nmaxmolatoms);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  widom_nmax = 0;
  ave_widom_chemical_potential = 0.0;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixWidomBiased::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix widom/biased command");

  // defaults

  exchmode = EXCHATOM;
  region_volume = 0;
  max_region_attempts = 1000;
  molecule_group = 0;
  molecule_group_bit = 0;
  molecule_group_inversebit = 0;
  exclusion_group = 0;
  exclusion_group_bit = 0;
  charge = 0.0;
  charge_flag = false;
  full_flag = false;
  energy_intra = 0.0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix widom mol", error);
      auto onemols = atom->get_molecule_by_id(arg[iarg+1]);
      if (onemols.size() == 0)
        error->all(FLERR,"Molecule template ID {} for fix widom does not exist", arg[iarg+1]);
      if (onemols.size() > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template {} for fix widom has multiple molecules; "
                       "will use only the first molecule", arg[iarg+1]);
      exchmode = EXCHMOL;
      onemol = onemols[0];
      iarg += 2;
    } 
      
      else if (strcmp(arg[iarg], "volumes") == 0) {
      this->read_volumes(arg[iarg]);
        
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      region = domain->get_region_by_id(arg[iarg+1]);
      if (!region) error->all(FLERR,"Region {} for fix widom does not exist",arg[iarg+1]);
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      charge = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"intra_energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      energy_intra = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix widom command");
  }
}

void FixWidomBiased::read_volumes(file){
}
