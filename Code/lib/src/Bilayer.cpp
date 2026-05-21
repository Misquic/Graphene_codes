#include "Bilayer.h"
#include "Constants.h" // freqquently and widely used variables used as parameters

//////////////////////////// local function declarations ///////////////////////

static inline double count1b(
  const double ni,
  const std::vector<double>& Capacities,
  const std::vector<double>& Voltages,
  const double nQ);

static inline double count1c(const double nC, const double nQ);

static inline double count1d(const std::vector<double>& Capacities);

inline double E_nl(double B, int abs_nl);

std::vector<double> prepare_E_nl(const double B, int L_max);

double prepare_part_sum(const std::vector<double>& E_nl_vec);

double count_E_0(
  const std::vector<double>& E_nl_vec,
  const double n0,
  const double B,
  const double part_sum);

// static double count_Vg(
//   const std::vector<double>& Cjg,
//   const std::vector<double>& Vj,
//   double E0,
//   const std::vector<double>& E_nl_vec,
//   double n0,
//   double B,
//   double part_sum);

static inline void count_n_V(double V, double Vg_o, double C, double Cg, double abs_ni, double sign_ni,
                      double &n, double &Vg);

///////////////////////////// method declarations //////////////////////////////

Bilayer::Bilayer(double nit, double nib,
                 double dt, double dg, double db,
                 double et, double eg, double eb):
  _nit{nit*Const::inv_cmsq2au}, _nib{nib*Const::inv_cmsq2au},
  _dt{std::abs(dt)*Const::nm2au}, _dg{std::abs(dg)*Const::nm2au}, _db{std::abs(db)*Const::nm2au},
  _et{std::abs(et)}, _eg{std::abs(eg)}, _eb{std::abs(eb)},
  _Ct{Const::eps_0*_et/_dt}, _Cg{Const::eps_0*_eg/_dg}, _Cb{Const::eps_0*_eb/_db}
{
};

Bilayer::Bilayer():
  _nit{Const::nit_default},
  _nib{Const::nib_default},
  _dt{Const::dt_default},
  _dg{Const::dg_default},
  _db{Const::db_default},
  _et{Const::et_default},
  _eg{Const::eg_default},
  _eb{Const::eb_default},
  _Ct{Const::Ct_default},
  _Cg{Const::Cg_default},
  _Cb{Const::Cb_default}
{
};

// pass Vt, Vb in eV
// return results in a.u.
std::pair<double, double> Bilayer::countDensitiesAndPotential(
  double Vt, double Vb, double* results) const
{
  LOG_VAR(Vt); //in eV
  LOG_VAR(Vb); // in eV

  constexpr double tol = 1e-4;
  Vt *= Const::eV2au; // Top gate voltage
  Vb *= Const::eV2au; // Bottom gate voltage

  double Vrange = Vt - Vb;

  double Vgt = Vb + 2.0/3.0*Vrange, Vgb = Vb + Vrange/3; // Graphene top/bottom voltage
  double Vgt_old = 0, Vgb_old = 0; // previous iteration Graphene top/bottom voltage
  double ngt      = 0, ngb      = 0; // Graphene top/bottom densities
  const double abs_nit = std::abs(_nit);
  const double abs_nib = std::abs(_nib);
  const double sign_nit = _nit >= 0. ? 1. : -1.;
  const double sign_nib = _nib >= 0. ? 1. : -1.;

  unsigned it = 0;
  bool converged = false;
  do
  {
    it++;
    // Vgt = Vgt_old;
    // Vgb = Vgb_old;
    // t
    count_n_V(Vt, Vgb_old, _Ct, _Cg, abs_nit, sign_nit, ngt, Vgt);
    // b
    count_n_V(Vb, Vgt_old, _Cb, _Cg, abs_nib, sign_nib, ngb, Vgb);


    if(it % 15 == 0)
    {
      if ((Vgt - Vgt_old)/Vgt_old < tol && (Vgb - Vgb_old)/Vgb_old < tol)
      {
          converged = true;
      }
      // printf("\r         \r%i", it);
    }
    Vgt_old = Vgt;
    Vgb_old = Vgb;
  } while ( !converged );

  // printf("iterations: %i\n", it);

  if(results)
  {
    results[0] = ngt / Const::inv_cmsq2au;
    results[1] = ngb / Const::inv_cmsq2au;
    results[2] = Vgt / Const::eV2au;
    results[3] = Vgb / Const::eV2au;
  }

  return {ngt, ngb};
};

// n0 is ni
// pass Vt, Vb, B in Si,
// return results in a.u.
// std::pair<double, double> Bilayer::countDensitiesAndPotential(
//   double Vt, double Vb, double B, double* results) const
// {
//   constexpr double tol = 1e-4;
//   LOG_VAR(Vt); //in eV
//   LOG_VAR(Vb); // in eV
//   LOG_VAR(B);  // in T
//   Vt *= Const::V2au;
//   Vb *= Const::V2au;
//   B  *= Const::T2au;

//   double Vrange = Vt - Vb;

//   // count E_nl -> array1D for const L const e const B various nl for sum later
//   LOG(std::vector<double> E_nl_vec = prepare_E_nl(B, Const::L_max));

//   LOG(double part_sum = prepare_part_sum(E_nl_vec));

//   // find E0 for n0/ni
//   // count basic n0 -> array2D for const L but variable B and variable E = x -> first use E_0 later use E_0 + eVg
//   LOG(double E0t = count_E_0(E_nl_vec, _nit, B, part_sum)); // quasi Fermi Energy top layer
//   LOG(double E0b = count_E_0(E_nl_vec, _nib, B, part_sum)); // quasi Fermi Energy bottom layer

//   std::cout << "E0t " << E0t << " E0b " << E0b << "\n";

//   // find Vg for E0
//   double Vgt = Vb + 2.0/3.0*Vrange, Vgb = Vb + Vrange/3; // Graphene top/bottom voltage
//   // double Vgt = 0, Vgb = 0;
//   double Vgt_old = 0, Vgb_old = 0;

//   // gates potentials for top flake
//   std::vector<double> Vjt = {Vt, Vgb}; // gate top, graphene bottom
//   // capacities for top flake
//   std::vector<double> Cjt = {_Ct, _Cg};

//   // gates potentials for bottom flake
//   std::vector<double> Vjb = {Vb, Vgt}; // gate bottom, graphene top
//   // capacities for bottom flake
//   std::vector<double> Cjb = {_Cb, _Cg};

//   uint it = 0;
//   const uint maxIt = 100;
//   bool converged = false;
//   do
//   {
//     dmsg("\nit " << it++ <<"\n");
//     // t
//     // set to value for other flake from prev iteration
//     Vjt[1] = Vgb_old;
//     LOG(Vgt = count_Vg(Cjt, Vjt, E0t, E_nl_vec, _nit, B, part_sum));

//     // b
//     // set to value for other flake from prev iteration
//     Vjb[1] = Vgt_old;
//     LOG(Vgb = count_Vg(Cjb, Vjb, E0b, E_nl_vec, _nib, B, part_sum));

//     double checkT = (Vgt - Vgt_old)/Vgt_old;
//     double checkB = (Vgb - Vgb_old)/Vgb_old;
//     LOG_VAR(checkT);
//     LOG_VAR(checkB);
//     if (checkT < tol && checkB < tol)
//     {
//       converged = true;
//     }
//     // printf("\r         \r%i", it);
//     Vgt_old = Vgt;
//     Vgb_old = Vgb;
//   } while ( !converged && it < maxIt);

//   if (it == maxIt)
//   {
//     std::cerr << "max nr of it, converged = " << (converged?"true":"false") << '\n';
//     exit(EXIT_FAILURE);
//   }
//   // obtained result for potentials at Graphene layers

//   // from eq.15 and
//   double ngt = _nit + _Cg*(Vgb - Vgt) + _Ct*(Vt - Vgt); // return this
//   double ngb = _nib + _Cg*(Vgt - Vgb) + _Cb*(Vb - Vgb); // return this

//   LOG_VAR(ngt / Const::inv_cmsq2au);
//   LOG_VAR(ngb / Const::inv_cmsq2au);
//   LOG_VAR(Vgt / Const::V2au);
//   LOG_VAR(Vgb / Const::V2au);
//   LOG_VAR(E0t / Const::eV2au);
//   LOG_VAR(E0b / Const::eV2au);

//   if(results)
//   {
//     results[0] = ngt / Const::inv_cmsq2au;
//     results[1] = ngb / Const::inv_cmsq2au;
//     results[2] = Vgt / Const::V2au;
//     results[3] = Vgb / Const::V2au;
//     results[4] = E0t / Const::eV2au;
//     results[5] = E0b / Const::eV2au;
//   }

//   return {ngt, ngb};

//   // TODO check if it works
// }

//////////////////////////// local function definitions ////////////////////////

static inline double count1b(
  const double ni,
  const std::vector<double>& Capacities,
  const std::vector<double>& Voltages,
  const double nQ)
{
  const double sum = (Capacities[0] * Voltages[0] +
                      Capacities[1] * Voltages[1]) /
                      Const::e;
  const double sqrtPart = signum(ni) * std::sqrt(2 * nQ * std::abs(ni));

  return ni + sum + sqrtPart;
}

static inline double count1c(const double nC, const double nQ)
{
  return signum(nC) * nQ * (1 - std::sqrt(1 + 2 * (std::abs(nC) / nQ)));
}

static inline double count1d(const std::vector<double>& Capacities)
{
  const double sum = (Capacities[0] + Capacities[1])/Const::e;

  return M_PI_2 * pow2(Const::hbar_v / Const::e * (sum));
}

// nl - number of Landau level nl = 0, -/+1, -/+2,...
inline double E_nl(double B, int abs_nl)
{
  return std::sqrt(2*B*Const::hbar_v*Const::hbar_v*abs_nl);
}

//precount E_nl for variable const L const e const B
std::vector<double> prepare_E_nl(const double B, int L_max)
{
  std::vector<double> E_nl_vec(std::abs(L_max));
  if(L_max > 0)
  {
    for(int i = 0; i < L_max; i++)
    {
      E_nl_vec[i] = E_nl(B, i);
    }
  }
  else
  {
    L_max = -L_max;
    for(int i = 0; i < L_max; i++)
    {
      E_nl_vec[i] = -E_nl(B, i);
    }
  }
  return E_nl_vec;
};

double prepare_part_sum(const std::vector<double>& E_nl_vec)
{
  double sum = 0;
  size_t size = E_nl_vec.size();
  for(size_t i = 0; i < size; i++)
  {
    sum += std::atan(E_nl_vec[i]*Const::inv_Lorentzian_par);
  }
  return sum;
};

// realises eq. nr 14
double count_E_0(
  const std::vector<double>& E_nl_vec,
  const double n0,
  const double B,
  const double part_sum)
{
  const double left    = n0*Const::h*M_PI/(4*B) - part_sum;
  const double E_low   = -1.0*Const::eV2au;                         //eV
  const double E_high  = 1.0*Const::eV2au;                          //eV
  dmsg("Searching E_0 in range: " << E_low << " - " << E_high << "\n");

  const size_t size = E_nl_vec.size();
  std::function<double(double)> func = [&E_nl_vec, size, left](double x)
    {
      double sum = 0;
      for(size_t i = 0; i < size; i++)
      {
        sum+= std::atan((x - E_nl_vec[i])*Const::inv_Lorentzian_par);
      }
      return sum - left;
    };

  return bisection(E_low, E_high, func);
};

// // realises eq. nr 16
// double count_Vg(double V_gate, double Vg_other, double C, double Cg, double E0,
//                 const std::vector<double>& E_nl_vec, double n0, double B, double part_sum)
// {
//   // left is the part without Vg
//   const double left = (n0 - part_sum*4*B/(Const::h*M_PI) + C*V_gate + Cg*Vg_other)*(Const::h*M_PI)/(4*B);
//   const double Vg_low = -15; //eV ?
//   const double Vg_high = 15; //eV ?
//   const double precalc = M_PI_4/B*(C + Cg);
//
//   const size_t size = E_nl_vec.size();
//   std::function<double(double)> func = [&E_nl_vec, size, left, E0, precalc](double Vg)
//     {
//       double sum = 0;
//       for(size_t i = 0; i < size; i++)
//       {
//         sum+= std::atan((E0 + Vg - E_nl_vec[i])*Const::inv_Lorentzian_par);
//         // sum += std::atan((Vg - E_nl_vec[i])*Const::inv_Lorentzian_par); // only when E_nl is precalculated
//       }
//       sum+= precalc*Vg;
//       return sum - left;
//     };
//
//   double Vg = bisection(Vg_low, Vg_high, func);
//   LOG_VAR(Vg_other)
//   LOG_VAR(Vg)
//   return Vg;
// };

// realises general eq. nr 16
// pass values in au units
// static double count_Vg(
//   const std::vector<double>& Cjg,
//   const std::vector<double>& Vj,
//   double E0,
//   const std::vector<double>& E_nl_vec,
//   double n0,
//   double B,
//   double part_sum)
// {
//   assert(Cjg.size() == Vj.size());
//   // partsum is sum_{nl}arctan(E_nl/eps)
//   // left is part without Vg
//   const double left = n0 - 4*B/(Const::h*M_PI) * part_sum;
//   // limits where to fing Vg
//   const double Vg_low = -15*Const::eV2au;
//   const double Vg_high = 15*Const::eV2au;

//   std::function<double(double)> bisectionFunction =
//     [&E_nl_vec, &Cjg, &Vj, left, E0, B](double Vg)
//       {
//         double sumE = 0, sumV = 0;

//         // sum for energies over landau levels
//         for(size_t nl = 0; nl < E_nl_vec.size(); nl++)
//         {
//           sumE += std::atan((E0 + Vg - E_nl_vec[nl])*Const::inv_Lorentzian_par);
//         }

//         // sum for potentials over "gates"
//         for (size_t j = 0; j < Cjg.size(); j++)
//         {
//           sumV += Cjg[j] * (Vj[j] - Vg);
//         }

//         return 4*B/(Const::h*M_PI)*sumE - sumV - left; // goes to 0 as Vg settles
//       };

//   double Vg = bisection(Vg_low, Vg_high, bisectionFunction);
//   LOG_VAR(Vg/Const::eV2au);
//   return Vg;
// }


inline void count_n_V(double V, double Vg_o, double C, double Cg, double abs_ni, double sign_ni,
                      double &n, double &Vg)
{
  double nQ      = M_PI_2 * pow2(Const::hbar_v * (C + Cg) );
  double sqrt_2_nQ_abs_ni = std::sqrt(2*nQ*abs_ni);
  double nC      = sign_ni * abs_ni + (C*V + Cg*Vg_o) + sign_ni * sqrt_2_nQ_abs_ni;
  double delta_n = signum(nC) * nQ * ( 1 - std::sqrt(1 + 2*std::abs(nC)/nQ));

  n = nC + delta_n;
  Vg = - (delta_n + sign_ni*sqrt_2_nQ_abs_ni) / (C + Cg);
}
