#ifndef FUNKCJE_H
#define FUNKCJE_H

#include "Utilities.h"

inline double signum(double x)
{
  return x > 0. ? 1. : -1.;
};

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
  std::pair<double, double> count_densities(double Vt, double Vb);

  /* returns nt, nb -> densities calculated with external magnetic field,
    provide Vt and Vb in eV and B in T */
  std::pair<double, double> count_densities(double Vt, double Vb, double B);

};

  std::vector<double>     prepare_E_nl(double B, int L_max);
  double prepare_part_sum(const std::vector<double>& E_nl_vec);
  double count_E_0(const std::vector<double>& E_nl_vec, double n0, double B, double partsum);

#endif


