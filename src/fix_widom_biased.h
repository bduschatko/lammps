/* -*- c++ -*- ----------------------------------------------------------

Widom insertion routine with sample biasing to system voids 

- Blake Duschatko

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(widom/biased,FixWidomBiased);
// clang-format on
#else

#ifndef LMP_FIX_WIDOM_BIASED_H
#define LMP_FIX_WIDOM_BIASED_H

#include "fix.h"
#include "fix_widom.h"

namespace LAMMPS_NS {

class FixWidomBiased : public FixWidom {
 public:

  FixWidomBiased(class LAMMPS *, int, char **);
  ~FixWidomBiased() override;

  // Adapted Methods

  void options(int narg, char **arg);
  double energy(int i, int itype, tagint imolecule, double *coord);

  // Pair Style ONLY Insertions
  //void attempt_atomic_insertion();
  void attempt_molecule_insertion();

  // FULL Energy Insertions 
  //void attempt_atomic_insertion_full();
  void attempt_molecule_insertion_full();

  // Custom Method

  void read_volumes(char *filename);

 private:

  // Array of volume centers and size 
  const Eigen::MatrixXd volumes;
  
};

}   
#endif
#endif
