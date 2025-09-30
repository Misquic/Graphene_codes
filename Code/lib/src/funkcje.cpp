#include "funkcje.h"
#include "Constants.h" // freqquently and widely used variables used as parameters

Bilayer::Bilayer(double nit, double nib,
                 double dt, double dg, double db,
                 double et, double eg, double eb):
  m_nit{nit*Const::inv_cmsq2au}, m_nib{nib*Const::inv_cmsq2au},
  m_dt{std::abs(dt)*Const::nm2au}, m_dg{std::abs(dg)*Const::nm2au}, m_db{std::abs(db)*Const::nm2au},
  m_et{std::abs(et)}, m_eg{std::abs(eg)}, m_eb{std::abs(eb)},
  m_Ct{Const::eps_0*m_et/m_dt}, m_Cg{Const::eps_0*m_eg/m_dg}, m_Cb{Const::eps_0*m_eb/m_db} 
{
};

inline void count_n_V(double V, double Vg_o, double C, double Cg, double abs_ni, double sign_ni,
                      double &n, double &Vg)
{
  double nQ      = M_PI_2 * pow2(Const::hv * (C + Cg) );
  double sqrt_2_nQ_abs_ni = std::sqrt(2*nQ*abs_ni);
  double nC      = sign_ni * abs_ni + (C*V + Cg*Vg_o) + sign_ni * sqrt_2_nQ_abs_ni;
  double delta_n = signum(nC) * nQ * ( 1 - std::sqrt(1 + 2*std::abs(nC)/nQ));
    
  n = nC + delta_n;
  Vg = - (delta_n + sign_ni*sqrt_2_nQ_abs_ni) / (C + Cg);
}

std::pair<double, double> Bilayer::count_densities(double Vt, double Vb)
{
  constexpr double tol = 1e-4;
  Vt *= Const::eV2au;
  Vb *= Const::eV2au;
    
  double Vgt     = 0, Vgb     = 0;
  double Vgt_old = 0, Vgb_old = 0;
  double nt      = 0, nb      = 0;
  const double abs_nit = std::abs(m_nit);
  const double abs_nib = std::abs(m_nib);
  const double sign_nit = m_nit >= 0. ? 1. : -1.;
  const double sign_nib = m_nib >= 0. ? 1. : -1.;

  unsigned it = 0;
  bool converged = false;
  do
  {
    it++;
    // Vgt = Vgt_old;
    // Vgb = Vgb_old;
    // t
    count_n_V(Vt, Vgb_old, m_Ct, m_Cg, abs_nit, sign_nit, nt, Vgt);
    // b
    count_n_V(Vb, Vgt_old, m_Cb, m_Cg, abs_nib, sign_nib, nb, Vgb);
    
    
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

  return {nt, nb};
};

// nl - number of Landau level nl = 0, -/+1, -/+2,...
inline double constexpr E_nl(double B, int abs_nl)
{
  return std::sqrt(2*B*Const::hv*Const::hv*abs_nl);
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
double count_E_0(const std::vector<double>& E_nl_vec,
                 const double n0,
                 const double B,
                 const double part_sum)
{
  const double left    = n0*Const::h*M_PI/(4*B) - part_sum;
  const double E_low   = -0.5*Const::eV2au;                         //eV
  const double E_high  = 0.5*Const::eV2au;                          //eV 
  dmsg("Searching in range: " << E_low << " - " << E_high << "\n");

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

// realises eq. nr 16
double count_Vg(double V, double Vg_o, double C, double Cg, double E0, const std::vector<double>& E_nl_vec, double n0, double B, double part_sum)
{
  const double left = (n0 - part_sum*4*B/(Const::h*M_PI) + C*V + Cg*Vg_o)*(Const::h*M_PI)/(4*B);
  const double Vg_low = -10; //eV ?
  const double Vg_high = 10; //eV ?
  const double precalc = M_PI_4/B*(C + Cg);

  const size_t size = E_nl_vec.size();
  std::function<double(double)> func = [&E_nl_vec, size, left, E0, precalc](double Vg)
    {
      double sum = 0;
      for(size_t i = 0; i < size; i++)
      {
        sum+= std::atan((E0 + Vg - E_nl_vec[i])*Const::inv_Lorentzian_par);
        // sum += std::atan((Vg - E_nl_vec[i])*Const::inv_Lorentzian_par); // only when E_nl is precalculated
      }
      sum+= precalc*Vg;
      return sum - left;
    };

  return bisection(Vg_low, Vg_high, func);
};

// n0 is ni
std::pair<double, double> Bilayer::count_densities(double Vt, double Vb, double B)
{
  constexpr double tol = 1e-4;
  Vt *= Const::eV2au;
  Vb *= Const::eV2au;
  B  *= Const::T2au;

  // count E_nl -> array1D for const L const e const B various nl for sum later
  std::vector<double> E_nl_vec = prepare_E_nl(B, Const::L_max);
  double part_sum = prepare_part_sum(E_nl_vec);

  // find E0 for n0/ni
  // count basic n0 -> array2D for const L but variable B and variable E = x -> first use E_0 later use E_0 + eVg 
  double E_0_t = count_E_0(E_nl_vec, m_nit, B, part_sum);
  double E_0_b = count_E_0(E_nl_vec, m_nib, B, part_sum);

  // find Vg for E0 
  double Vgt     = 0, Vgb     = 0;
  double Vgt_old = 0, Vgb_old = 0;

  unsigned it = 0;
  bool converged = false;
  do
  {
    it++;
    // t
    Vgt = count_Vg(Vt, Vgb_old, m_Ct, m_Cg, E_0_t, E_nl_vec, m_nit, B, part_sum);
    // b
    Vgb = count_Vg(Vb, Vgt_old, m_Cb, m_Cg, E_0_b, E_nl_vec, m_nib, B, part_sum);

    if ((Vgt - Vgt_old)/Vgt_old < tol && (Vgb - Vgb_old)/Vgb_old < tol)
    {
      converged = true;
    }
    // printf("\r         \r%i", it);
    Vgt_old = Vgt;
    Vgb_old = Vgb;
  } while ( !converged );

  // TODO return n!
  // TODO check if it works
  // TODO get results
}

//
