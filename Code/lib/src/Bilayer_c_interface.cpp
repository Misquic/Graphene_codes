#include "Bilayer.h"
#include "Constants.h"

Bilayer* Bilayer_constructor(double nit, double nib,
                             double dt, double dg, double db,
                             double et, double eg, double eb)
{
  return new Bilayer(nit, nib, dt, dg, db, et, eg, eb);
};

void Bilayer_destructor(Bilayer* this_ptr)
{
  delete this_ptr;
  this_ptr = NULL;
};

void Bilayer_countDensities(const Bilayer* const this_ptr,
                            double Vt, double Vb,
                            double* const nt, double* const nb)
{
  auto[nt_temp, nb_temp] = this_ptr->countDensities(Vt, Vb);
  *nt = nt_temp;
  *nb = nb_temp;
};

void Bilayer_countDensities_B(const Bilayer* const this_ptr,
                              double Vt, double Vb, double B,
                              double* const nt, double* const nb)
{
  auto[nt_temp, nb_temp] = this_ptr->countDensities(Vt, Vb, B);
  *nt = nt_temp;
  *nb = nb_temp;
};

double count_E_0(double n0, double B){
  std::vector<double> E_nl_vec(std::abs(Const::L_max));
  
  if(n0 > 0)
  {
    E_nl_vec = prepare_E_nl(B, std::abs(Const::L_max));
  }
  else
  {
    E_nl_vec = prepare_E_nl(B, -std::abs(Const::L_max));
  }

  return count_E_0(E_nl_vec, n0, B, prepare_part_sum(E_nl_vec))/Const::eV2au;
};
