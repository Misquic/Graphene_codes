#include "Bilayer.h"
#include "BilayerPositionGeneration.h"
#include "Constants.h"

Bilayer* Bilayer_constructor(double nit, double nib,
                             double dt, double dg, double db,
                             double et, double eg, double eb)
{
  return new Bilayer(nit, nib, dt, dg, db, et, eg, eb);
};

Bilayer* Bilayer_constructor_default()
{
  Bilayer* bilayer_p = new Bilayer();
  dmsg("bilayer_p " << bilayer_p << "\n");
  return bilayer_p;
};

void Bilayer_destructor(Bilayer* this_ptr)
{
  delete this_ptr;
  this_ptr = NULL;
};

void Bilayer_countDensities(
  const Bilayer* const this_ptr,
  double Vt,
  double Vb,
  double* const nt_p,
  double* const nb_p)
{
  auto[nt_temp, nb_temp] = this_ptr->countDensitiesAndPotential(Vt, Vb);
  *nt_p = nt_temp;
  *nb_p = nb_temp;
};

void Bilayer_countEnergiesAndPotential(
  const Bilayer* const this_ptr,
  const double Vt,
  const double Vb,
  double* const Vgt_p,
  double* const Vgb_p)
{
  double results[4] = {};

  this_ptr->countDensitiesAndPotential(Vt, Vb, results);

  *Vgt_p = results[2];
  *Vgb_p = results[3];
}

void Bilayer_countDensities_B(
  const Bilayer* const this_ptr,
  double Vt,
  double Vb,
  double B,
  double* const nt_p,
  double* const nb_p)
{
  Bilayer::resultsB results = this_ptr->countDensitiesAndPotential(Vt, Vb, B);
  *nt_p = results.nt;
  *nb_p = results.nb;
};

void Bilayer_countEnergiesAndPotential_B(
  const Bilayer* const this_ptr,
  const double Vt,
  const double Vb,
  const double B,
  double* const Vgt_p,
  double* const Vgb_p,
  double* const E0t_p,
  double* const E0b_p)
{
  Bilayer::resultsB results = this_ptr->countDensitiesAndPotential(Vt, Vb, B);

  *Vgt_p = results.Vgt;
  *Vgb_p = results.Vgb;
  *E0t_p = results.E0t;
  *E0b_p = results.E0b;
}

void Bilayer_countAll_B(
  const Bilayer* const this_ptr,
  const double Vt,
  const double Vb,
  const double B,
  double* const Vgt_p,
  double* const Vgb_p,
  double* const E0t_p,
  double* const E0b_p,
  double* const nt_p,
  double* const nb_p)
{
  Bilayer::resultsB results = this_ptr->countDensitiesAndPotential(Vt, Vb, B);

  *Vgt_p = results.Vgt;
  *Vgb_p = results.Vgb;
  *E0t_p = results.E0t;
  *E0b_p = results.E0b;
  *nt_p = results.nt;
  *nb_p = results.nb;
}

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

extern "C" void Bilayer_generatePositions3D(
  const int nX,
  const int nY,
  const double foldRadius,
  const int cutLead,
  const int leadWidth,
  const double** const data_pp,
  int* const size_p)
{
  ParamsS params =
  {
    .nX = (size_t)nX,
    .nY = (size_t)nY,
    .foldRadius = foldRadius,
    .cutLead = (size_t)cutLead,
    .leadWidth = (size_t)leadWidth
  };

  makeParamsRight(params);

  static std::vector<Vec3S> positions = generatePositions(params);

#ifdef DEBUG
  savePositions(params.saveFileName, positions);
#endif

  *data_pp = &positions[0].x;
  *size_p = positions.size();
}
