#include "Bilayer.h"
#include "Constants.h" // frequently and widely used variables used as parameters
#include <array>
#include "Utilities.h"

extern Rnd rnd;
//////////////////////////// local function declarations ///////////////////////

static inline std::array<double, 2> prepareTop(
  const double topValue, const double grapheneValue);

static inline std::array<double, 2> prepareBot(
  const double grapheneValue, const double bottomValue);

static inline size_t getLevelIndex(const int level);

static inline std::array<double, Const::Nl>
  preparelandauLevelsEnergiesTimesInvLorenzianPar(const double Bz);

static inline double countN(
  const double n0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double Vg);

static inline double landauLevelEnergy(const int level, const double Bz);

static inline double count4eBz_hPi(const double Bz);

static double countE0(
  const double Bz,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar);

static double countVg(
  const double Bz,
  const double E0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar);

static void checkE0(const double E0,
                    const double B_au,
                    const double _ni,
                    const std::array<double,
                    Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar);

///////////////////////////// method declarations //////////////////////////////

/* returns {nt, nb, Vgt, Vgb, E0t, E0b} -> densities calculated with external magnetic field,
    provide Vt and Vb in eV and B in T */
// TODO compare equations with original paper not tdmac5536supp1.pdf
Bilayer::resultsB Bilayer::countDensitiesAndPotential(
  double Vt, double Vb, double B) const
{
  const double Vt_au = Vt * Const::V2au;
  const double Vb_au = Vb * Const::V2au;
  const double B_au = B * Const::T2au;

  // initial values and conversion to au
  double Vgt = clamp((Vb_au + 2./3. * (Vt_au - Vb_au)) * 0.01,
                     Const::VgMin,
                     Const::VgMax);
  double Vgb = clamp((Vb_au + 1./3. * (Vt_au - Vb_au)) * 0.01,
                     Const::VgMin,
                     Const::VgMax);

  // precalculate some values
  // std::cout << "precalc\n";
  const std::array<double, Const::Nl> landauLevelsEnergiesTimesInvLorenzianPar =
    preparelandauLevelsEnergiesTimesInvLorenzianPar(B_au);

  const double E0t = countE0(B_au, _nit, landauLevelsEnergiesTimesInvLorenzianPar);
  // const double E0t = -0.0319331884 * Const::eV2au;
  const double E0b = countE0(B_au, _nib, landauLevelsEnergiesTimesInvLorenzianPar);
  // const double E0b = -0.0319331884 * Const::eV2au;

  bool converged = false;
  uint8_t it = 0;

  do
  {
    it++;
    auto[newVgt, newVgb, newConverged] =
      iterationB(
        Vt_au, Vb_au, Vgt, Vgb, B_au, E0t, E0b, landauLevelsEnergiesTimesInvLorenzianPar);

    Vgt = Const::VgAlpha * newVgt + (1. - Const::VgAlpha) * Vgt;
    Vgb = Const::VgAlpha * newVgb + (1. - Const::VgAlpha) * Vgb;
    converged = newConverged;

  } while (!converged && it < Const::maxIterationsB);

  assertVerbose(converged, "IterationB didn't converge, Vgt=" << Vgt << " Vgb=" << Vgb << '\n');

  // count densities for final voltages
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);
  const std::array<double, 2> voltagesT = prepareTop(Vt_au, Vgb);
  const std::array<double, 2> voltagesB = prepareBot(Vgt, Vb_au);

  const double nt = countN(_nit, capacitiesT, voltagesT, Vgt);
  const double nb = countN(_nib, capacitiesB, voltagesB, Vgb);

#ifdef DEBUG
  checkE0(E0t, B_au, _nit, landauLevelsEnergiesTimesInvLorenzianPar);
  checkE0(E0b, B_au, _nib, landauLevelsEnergiesTimesInvLorenzianPar);
#endif

  resultsB results = {
    .nt = nt / Const::inv_cmsq2au,
    .nb = nb / Const::inv_cmsq2au,
    .Vgt = Vgt / Const::V2au,
    .Vgb = Vgb / Const::V2au,
    .E0t = E0t / Const::eV2au,
    .E0b = E0b / Const::eV2au
  };

  // std::printf("-E0t - eVgt = %f \t -E0b - eVgb = %f\n",
  //            -results.E0t - results.Vgt,
  //            -results.E0b - results.Vgb);

  return results;
};

// returns Vgt, Vgb, converged
std::tuple<double, double, bool> Bilayer::iterationB(
  const double Vt,
  const double Vb,
  const double prevVgt,
  const double prevVgb,
  const double Bz,
  const double E0t,
  const double E0b,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar) const
{
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);

  std::array<double, 2> voltagesT = prepareTop(Vt, prevVgb);
  std::array<double, 2> voltagesB = prepareBot(prevVgt, Vb);

  double Vgt = 0;
  double Vgb = 0;

  Vgt =
    countVg(
      Bz,
      E0t,
      capacitiesT,
      voltagesT,
      _nit,
      landauLevelsEnergiesTimesInvLorenzianPar);

  voltagesB[0] = Vgt;
  Vgb =
    countVg(
      Bz,
      E0b,
      capacitiesB,
      voltagesB,
      _nib,
      landauLevelsEnergiesTimesInvLorenzianPar);

  const double absoluteErrorT = std::abs(Vgt - prevVgt);
  const double absoluteErrorB = std::abs(Vgb - prevVgb);

  const bool convergedT =
    absoluteErrorT < Const::VgAbsTol + Const::VgRelTol * std::abs(Vgt);
  const bool convergedB =
    absoluteErrorB < Const::VgAbsTol + Const::VgRelTol * std::abs(Vgb);

  const bool converged = convergedT && convergedB;

  return {Vgt, Vgb, converged};
}

//////////////////////////// local function definitions ////////////////////////

// above top layer, below top layer
static inline std::array<double, 2> prepareTop(
  const double topValue, const double grapheneValue)
{
  return {topValue, grapheneValue};
}

// above bottom layer, below bottom layer
static inline std::array<double, 2> prepareBot(
  const double grapheneValue, const double bottomValue)
{
  return {grapheneValue, bottomValue};
}

static inline size_t getLevelIndex(const int level)
{
  assert(level >= -Const::L_max && level <= Const::L_max);
  return Const::L_max + level;
}

static inline std::array<double, Const::Nl>
  preparelandauLevelsEnergiesTimesInvLorenzianPar(const double Bz)
{
  std::array<double, Const::Nl> landauLevelsEnergiesTimesInvLorenzianPar{0};

  for (int level = -Const::L_max; level <= Const::L_max; level++)
  {
    // dmsg("index: " << getLevelIndex(level) << " landauLevelEnergy: " << landauLevelEnergy(level, Bz) << '\n');
    landauLevelsEnergiesTimesInvLorenzianPar[getLevelIndex(level)] =
      landauLevelEnergy(level, Bz) * Const::inv_Lorentzian_par;
  }

  return landauLevelsEnergiesTimesInvLorenzianPar;
}

static inline double countN(
  const double n0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double Vg)
{
  const double nG = (capacities[0] * (voltages[0] - Vg) +
                     capacities[1] * (voltages[1] - Vg)) / Const::e;

  return n0 + nG;
}

// level = nl in paper
static inline double landauLevelEnergy(const int level, const double Bz)
{
  return signum(level) *
    std::sqrt(
      2 * Const::e * Bz * Const::hbar * Const::v * Const::v * std::abs(level));
}

static inline double count4eBz_hPi(const double Bz)
{
  return (4. * Const::e * Bz) / (Const::h * M_PI);
}

static double countE0(
  const double Bz,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, expr, &landauLevelsEnergiesTimesInvLorenzianPar]
  (const double E0)
  {
    const double E0timesInvLorenzianPar = E0 * Const::inv_Lorentzian_par;

    // double sum = std::atan(landauLevelsEnergiesTimesInvLorenzianPar[0] + 0); // 0th level
    double sum = 0;
    for (int level = -Const::L_max; level <= Const::L_max; level++)
    {
      double EnlTimesInvLorPar = landauLevelsEnergiesTimesInvLorenzianPar[getLevelIndex(level)];
      sum += (std::atan(E0timesInvLorenzianPar - EnlTimesInvLorPar));
      //  +
              // std::atan(E0timesInvLorenzianPar + EnlTimesInvLorPar));
      // E_{nl} = -E_{-nl}
    }

    // sum += sumAtan_EnlTimesInvLorenzianPar; // sumAtan = 0
    return expr * sum - n0;
  };

  // std::cout << "E0\n";
  return bisection(Const::E0Min, Const::E0Max, func);
}

static double countVg(
  const double Bz,
  const double E0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, &capacities, &voltages, expr, Bz, E0,
    &landauLevelsEnergiesTimesInvLorenzianPar]
  (double Vg)
  {
    const double left = countN(n0, capacities, voltages, Vg);
    const double shiftedE = (E0 + Const::e * Vg) * Const::inv_Lorentzian_par;

    // double right = std::atan(landauLevelsEnergiesTimesInvLorenzianPar[0] + 0); // 0th level
    double right = 0;
    for (int level = -Const::L_max; level <= Const::L_max; level++)
    {
      double EnlTimesInvLorPar = landauLevelsEnergiesTimesInvLorenzianPar[getLevelIndex(level)];
      right += (std::atan(shiftedE - EnlTimesInvLorPar));
      //  +
                // std::atan(shiftedE + EnlTimesInvLorPar));
      // E_{nl} = -E_{-nl}
    }

    // right += sumAtan_EnlTimesInvLorenzianPar // sumAtan = 0;

    return expr*right - left;
  };

  // std::cout << "Vg\n";
  return bisection(Const::VgMin, Const::VgMax, func);
}

static void checkE0(const double E0,
                    const double B_au,
                    const double _ni,
                    const std::array<double,
                    Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar)
{
  double sum = 0;
  const double expr = count4eBz_hPi(B_au);
  (void)landauLevelsEnergiesTimesInvLorenzianPar;
  // const double E0timesInvLorenzianPar = E0 * Const::inv_Lorentzian_par;
  // for (double EnlTimesInvLorenzianPar: landauLevelsEnergiesTimesInvLorenzianPar)
  // {
  //   sum += std::atan(E0timesInvLorenzianPar - EnlTimesInvLorenzianPar);
  // }

  // sum += sumAtan_EnlTimesInvLorenzianPar;

  for (int level = -Const::L_max; level <= Const::L_max; level++)
  {
    sum += std::atan((E0 - landauLevelEnergy(level, B_au)) * Const::inv_Lorentzian_par) + 0;
          //  std::atan(landauLevelEnergy(level, B_au) * Const::inv_Lorentzian_par);
  }

  double check = expr * sum;

  std::printf("_ni = %e, check = %e, diff = %e, rel = %f%% \n", _ni, check, std::abs(_ni - check), std::abs(_ni - check)/check * 100);
}

