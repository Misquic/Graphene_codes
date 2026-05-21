#ifndef BILAYER_H
#define BILAYER_H

#include "Utilities.h"
#include "Constants.h"

class Bilayer
{
private:
  const double _nit, _nib; // initial densities n0
  const double _dt, _dg, _db;
  const double _et, _eg, _eb;

  const double _Ct, _Cg, _Cb;

  std::tuple<double, double, bool> iterationB(
    const double Vt,
    const double Vb,
    const double prevVgt,
    const double prevVgb,
    const double Bz,
    const double E0t,
    const double E0b,
    const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar,
    const double sumAtan) const;

public:
  // provide densities <n> in 1/cm^2, thicknesses <d> in nm, permitivites <e> in -
  Bilayer(double nit, double nib,
          double dt, double dg, double db,
          double et, double eg, double eb);

  Bilayer();

  Bilayer(Bilayer&& other) = delete;
  Bilayer(Bilayer& other) = delete;
  Bilayer operator=(Bilayer&& other) = delete;
  Bilayer operator=(Bilayer& other) = delete;

  /* returns nt, nb -> densities calculated without external magnetic field,
     returns [nt, nb, Vgt, Vgb] in results* if given
     provide Vt and Vb in eV */
  std::pair<double, double> countDensitiesAndPotential(
    double Vt, double Vb, double* results = nullptr) const;

  typedef struct resultsB
  {
    double nt;
    double nb;
    double Vgt;
    double Vgb;
    double E0t;
    double E0b;
  } resultsB;

  // TODO compare equations with original paper not tdmac5536supp1.pdf
  /* returns {nt, nb, Vgt, Vgb, E0t, E0b} -> densities calculated with external magnetic field,
      provide Vt and Vb in eV and B in T */
  resultsB countDensitiesAndPotential(
    double Vt, double Vb, double B) const;

};

  // visible outside for testing
  std::vector<double> prepare_E_nl(double B, int L_max);
  double prepare_part_sum(const std::vector<double>& E_nl_vec);
  double count_E_0(const std::vector<double>& E_nl_vec, double n0, double B, double partsum);

/////////////////// c interface ///////////////////
extern "C"
{
  Bilayer* Bilayer_constructor(double nit, double nib,
                               double dt, double dg, double db,
                               double et, double eg, double eb);
  Bilayer* Bilayer_constructor_default();

  void Bilayer_destructor(Bilayer* this_ptr);

  void Bilayer_countDensities(
    const Bilayer* const this_ptr,
    double Vt,
    double Vb,
    double* const nt_p,
    double* const nb_p);

  void Bilayer_countEnergiesAndPotential(
    const Bilayer* const this_ptr,
    const double Vt,
    const double Vb,
    double* const Vgt_p,
    double* const Vgb_p);

  void Bilayer_countDensities_B(
    const Bilayer* const this_ptr,
    double Vt,
    double Vb,
    double B,
    double* const nt_p,
    double* const nb_p);

  void Bilayer_countEnergiesAndPotential_B(
    const Bilayer* const this_ptr,
    const double Vt,
    const double Vb,
    const double B,
    double* const Vgt_p,
    double* const Vgb_p,
    double* const E0t_p,
    double* const E0b_p);

  void Bilayer_countAll_B(
    const Bilayer* const this_ptr,
    const double Vt,
    const double Vb,
    const double B,
    double* const Vgt_p,
    double* const Vgb_p,
    double* const E0t_p,
    double* const E0b_p,
    double* const nit_p,
    double* const nib_p);

  double count_E_0(double n0, double B);
}

#endif
