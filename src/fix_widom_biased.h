/* -*- c++ -*- ----------------------------------------------------------

Parent members that require no changes

 public:
  FixWidom(class LAMMPS *, int, char **);
  ~FixWidom() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;

  void attempt_atomic_insertion();
  void attempt_molecule_insertion();

  void attempt_atomic_insertion_full();
  void attempt_molecule_insertion_full();
  double energy(int, int, tagint, double *);
  double molecule_energy(tagint);
  double energy_full();
  double compute_vector(int) override;
  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void grow_molecule_arrays(int);

 protected:
  int molecule_group, molecule_group_bit;
  int molecule_group_inversebit;
  int exclusion_group, exclusion_group_bit;
  int nwidom_type, nevery, seed;
  int ninsertions;
  int exchmode;            // exchange ATOM or MOLECULE
  class Region *region;    // widom region
  char *idregion;          // widom region id
  bool charge_flag;        // true if user specified atomic charge
  bool full_flag;          // true if doing full system energy calculations

  int natoms_per_molecule;    // number of atoms in each inserted molecule
  int nmaxmolatoms;           // number of atoms allocated for molecule arrays

  double ave_widom_chemical_potential;

  int widom_nmax;
  int max_region_attempts;
  double gas_mass;
  double insertion_temperature;
  double beta, volume;
  double charge;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double region_xlo, region_xhi, region_ylo, region_yhi, region_zlo, region_zhi;
  double region_volume;
  double energy_stored;    // full energy of old/current configuration
  double *sublo, *subhi;
  double **cutsq;
  double **molcoords;
  double *molq;
  imageint *molimage;
  imageint imagezero;

  double energy_intra;

  class Pair *pair;

  class RanPark *random_equal;

  class Atom *model_atom;

  class Molecule *onemol;
  int triclinic;    // 0 = orthog box, 1 = triclinic

  class Compute *c_pe;

  void options(int, char **);

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(widom/biased,FixWidomBiased);
// clang-format on
#else

#ifndef LMP_FIX_WIDOM_BIASED_H
#define LMP_FIX_WIDOM_BIASED_H

#include "fix_widom.h"


namespace LAMMPS_NS {

class FixWidomBiased : public FixWidom {

 public:

  FixWidomBiased(class LAMMPS *, int, char **); // constructor
  ~FixWidomBiased() override; // destructor 

  // Parent function must get called at beginning of Verlet run 
  // Only thing we need is to set the full flag for allegro 
  // It must be overriden by default but can call the parent from within 
  void init() override;

  // This calls the insertion routines and will call the parent ones if
  // not completely overriden by this method 
  void pre_exchange() override;

  // Local versions, not overriden but called locally
  void attempt_atomic_insertion();
  void attempt_molecule_insertion();

  void attempt_atomic_insertion_full();
  void attempt_molecule_insertion_full();

  void attempt_molecule_insertion_allegro();

 protected:

  // Does NOT override Fix Widom
  void options(int, char **);

  void read_volumes(char *filename);

  // New members
  bool allegro_flag = false;
  bool dynamic_volume;
  int nvols;
  int nphi = 1;
  int max_inserts = 1000000;
  double spacing = 1.0;
  double solvent_buffer = 1.0;
  double *volume_buffer;

  std::vector<std::vector<double>> volumes;

};

}    // namespace LAMMPS_NS

#endif
#endif
