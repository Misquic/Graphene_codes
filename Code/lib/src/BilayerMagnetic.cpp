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
  preparelandauLevelsEnergies(const double Bz);

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
                    const double _ni);

///////////////////////////// method declarations //////////////////////////////

/* returns {nt, nb, Vgt, Vgb, E0t, E0b} -> densities calculated with external magnetic field,
    provide Vt and Vb in eV and B in T */
Bilayer::resultsB Bilayer::countDensitiesAndPotential(
  double Vt, double Vb, double B) const
{
  const double Vt_au = Vt * Const::V2au;
  const double Vb_au = Vb * Const::V2au;
  const double B_au = B * Const::T2au;

  double Vgt = 0;
  double Vgb = 0;

  // precalculate some values
  const std::array<double, Const::Nl> landauLevelsEnergies =
    preparelandauLevelsEnergies(B_au);
  const double E0t = countE0(B_au, _nit, landauLevelsEnergies);
  const double E0b = countE0(B_au, _nib, landauLevelsEnergies);

  bool converged = false;
  uint16_t it = 0;

  do
  {
    it++;
    auto[newVgt, newVgb, newConverged] =
      iterationB(Vt_au, Vb_au, Vgt, Vgb, B_au, E0t, E0b, landauLevelsEnergies);

    Vgt = Const::VgAlpha * newVgt + (1. - Const::VgAlpha) * Vgt;
    Vgb = Const::VgAlpha * newVgb + (1. - Const::VgAlpha) * Vgb;
    converged = newConverged;

  } while (!converged && it < Const::maxIterationsB);

  dmsg(it);
  assertVerbose(converged, "IterationB didn't converge, Vgt=" << Vgt << " Vgb=" << Vgb << '\n');

  // count densities for final voltages
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);
  const std::array<double, 2> voltagesT = prepareTop(Vt_au, Vgb);
  const std::array<double, 2> voltagesB = prepareBot(Vgt, Vb_au);

  const double nt = countN(_nit, capacitiesT, voltagesT, Vgt);
  const double nb = countN(_nib, capacitiesB, voltagesB, Vgb);

  resultsB results = {
    .nt = nt / Const::inv_cmsq2au,
    .nb = nb / Const::inv_cmsq2au,
    .Vgt = Vgt / Const::V2au,
    .Vgb = Vgb / Const::V2au,
    .E0t = E0t / Const::eV2au,
    .E0b = E0b / Const::eV2au
  };

  return results;
};

double checkDiffVg(const double Vg,
                   const double Bz,
                   const double E0,
                   const std::array<double, 2>& capacities,
                   const std::array<double, 2>& voltages,
                   const double n0,
                   const std::array<double, Const::Nl>& landauLevelsEnergies)
{
  const double expr = count4eBz_hPi(Bz);

  const double left = countN(n0, capacities, voltages, Vg);

  const double shiftedE = E0 + Const::e * Vg;
  double right = std::atan((shiftedE - landauLevelsEnergies[0]) * Const::inv_Lorentzian_par); // 0th level
  for (int level = 1; level <= Const::L_max; level++)
  {
    double Enl = landauLevelsEnergies[getLevelIndex(level)];
    right += (std::atan((shiftedE - Enl) * Const::inv_Lorentzian_par) +
              std::atan((shiftedE + Enl) * Const::inv_Lorentzian_par));
  }

  return left - right * expr;
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
  const std::array<double, Const::Nl>& landauLevelsEnergies) const
{
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);

  std::array<double, 2> voltagesT = prepareTop(Vt, prevVgb);
  std::array<double, 2> voltagesB = prepareBot(prevVgt, Vb);

  double Vgt = 0;
  double Vgb = 0;

  Vgb = countVg(Bz,
                E0b,
                capacitiesB,
                voltagesB,
                _nib,
                landauLevelsEnergies);

  voltagesT[1] = Vgb;

  Vgt = countVg(Bz,
                E0t,
                capacitiesT,
                voltagesT,
                _nit,
                landauLevelsEnergies);

  voltagesB[0] = Vgt;

  const double absoluteErrorT = std::abs(checkDiffVg(Vgt, Bz, E0t, capacitiesT, voltagesT, _nit, landauLevelsEnergies));
  const double absoluteErrorB = std::abs(checkDiffVg(Vgb, Bz, E0b, capacitiesB, voltagesB, _nib, landauLevelsEnergies));

  dmsg(Vgt / Const::eV2au << " " << absoluteErrorT << " " << Vgb / Const::eV2au << " " << absoluteErrorB);
  const bool convergedT =
    absoluteErrorT < Const::VgAbsTol;
  const bool convergedB =
    absoluteErrorB < Const::VgAbsTol;

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
  return level;
}

static inline std::array<double, Const::Nl>
  preparelandauLevelsEnergies(const double Bz)
{
  std::array<double, Const::Nl> landauLevelsEnergies{0};

  for (int level = 0; level <= Const::L_max; level++)
  {
    landauLevelsEnergies[getLevelIndex(level)] = landauLevelEnergy(level, Bz);
  }

  return landauLevelsEnergies;
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
  const std::array<double, Const::Nl>& landauLevelsEnergies)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, expr, &landauLevelsEnergies]
  (const double E0)
  {
    const double E0timesInvLorenzianPar = E0 * Const::inv_Lorentzian_par;

    double sum = std::atan((E0 - landauLevelsEnergies[0]) * Const::inv_Lorentzian_par); // 0th level
    for (int level = 1; level <= Const::L_max; level++)
    {
      const double Enl = landauLevelsEnergies[getLevelIndex(level)];
      sum += (std::atan((E0 - Enl) * Const::inv_Lorentzian_par) +
              std::atan((E0 + Enl) * Const::inv_Lorentzian_par)); // E_{nl} = -E_{-nl}
    }

    return expr * sum - n0;
  };

  return bisection(Const::E0Min, Const::E0Max, func);
}

static double countVg(
  const double Bz,
  const double E0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergies)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, &capacities, &voltages, expr, Bz, E0, &landauLevelsEnergies]
  (double Vg)
  {
    const double left = countN(n0, capacities, voltages, Vg);
    const double shiftedE = E0 + Const::e * Vg;

    double right = std::atan((shiftedE - landauLevelsEnergies[0]) * Const::inv_Lorentzian_par); // 0th level
    for (int level = 1; level <= Const::L_max; level++)
    {
      double Enl = landauLevelsEnergies[getLevelIndex(level)];
      right += (std::atan((shiftedE - Enl) * Const::inv_Lorentzian_par) +
                std::atan((shiftedE + Enl) * Const::inv_Lorentzian_par)); // E_{nl} = -E_{-nl}
    }

    return expr*right - left;
  };

  return bisection(Const::VgMin, Const::VgMax, func);
}

static void checkE0(const double E0,
                    const double B_au,
                    const double _ni)
{
  double sum = 0;
  const double expr = count4eBz_hPi(B_au);

  for (int level = -Const::L_max; level <= Const::L_max; level++)
  {
    sum += std::atan((E0 - landauLevelEnergy(level, B_au)) * Const::inv_Lorentzian_par) + 0;
  }

  double check = expr * sum;

  std::printf("_ni = %e, check = %e, diff = %e, rel = %f%% \n", _ni, check, std::abs(_ni - check), std::abs(_ni - check)/check * 100);
}

