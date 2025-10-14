#ifndef BILAYER_H
#define BILAYER_H

#include "Utilities.h"

class Bilayer
{
private:
  const double m_nit, m_nib;
  const double m_dt, m_dg, m_db;
  const double m_et, m_eg, m_eb;

  const double m_Ct, m_Cg, m_Cb;

public: 
  // provide densities <n> in 1/cm^2, thicknesses <d> in nm, permitivites <e> in -
  Bilayer(double nit, double nib,
          double dt, double dg, double db,
          double et, double eg, double eb);

  Bilayer(Bilayer&& other) = delete;
  Bilayer(Bilayer& other) = delete;
  Bilayer operator=(Bilayer&& other) = delete;
  Bilayer operator=(Bilayer& other) = delete;

  /* returns nt, nb -> densities calculated without external magnetic field, 
    provide Vt and Vb in eV */
  std::pair<double, double> countDensities(double Vt, double Vb) const;

  /* returns nt, nb -> densities calculated with external magnetic field,
    provide Vt and Vb in eV and B in T */
  std::pair<double, double> countDensities(double Vt, double Vb, double B) const;

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
  void Bilayer_destructor(Bilayer* this_ptr);
  void Bilayer_countDensities(const Bilayer* const this_ptr,
                              double Vt, double Vb,
                              double* const nt, double* const nb);
  void Bilayer_countDensities_B(const Bilayer* const this_ptr,
                                double Vt, double Vb, double B,
                                double* const nt, double* const nb);

  double count_E_0(double n0, double B);
}

#endif
