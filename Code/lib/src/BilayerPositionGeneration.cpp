#include "BilayerPositionGeneration.h"

//////////////// local function declarations ////////////////

static Vec3S genFlatPos(const size_t i, const size_t j, const uint8_t atom);

static Vec3S genUpperPos(
  const size_t j,
  const Vec3S& prevPos,
  const double foldRadius,
  const bool firstInUpper);

static Vec3S genFold(
  const size_t i,
  const size_t j,
  const std::vector<Vec3S>& positions,
  const double foldRadius,
  const ParamsS& params);

//////////////// global function definitions ////////////////

std::vector<Vec3S>& generatePositions(const ParamsS& params)
{
  static std::vector<Vec3S> positions;

  const double R = params.foldRadius;
  const double rowSpacing = Const::vecsArmchair[1][1];
  const size_t nInFold = static_cast<size_t>(std::round(M_PI * R / rowSpacing));
  positions.clear();
  positions.reserve(params.nX * (params.nY * 2 + nInFold) * 2);

  double maxZ = 0;
  size_t lastY = params.nY * 2 + nInFold;

  for (uint8_t atom = 0; atom < 2; atom++)
  {
    for (size_t i = 0; i < params.nX; i++)
    {
      Vec3S lastFoldPos = {};
      bool hasFoldPos = false;
      Vec3S prevUpperPos = {};
      bool hasUpperPos = false;

      for (size_t j = 0; j < lastY; j++)
      {
        bool add = true;
        // cut bad edges
        if ((j == 0 && atom == 0) || (j == lastY - 1 && atom == 1))
        {
          add = false;
        }

        // cut bad corners
        if ((i <= 2 && j == 0 && atom == 1) ||
            (i <= 1 && j <= 1 && atom == 0) ||
            ((i >= params.nX - 3) && (j == (lastY - 1)) && (atom == 0)) ||
            ((i >= params.nX - 2) && (j == (lastY - 2)) && (atom == 1)))
        {
          add = false;
        }

        // cut for X lower lead
        if ((i <= 1) &&
            ((j <= (params.nY - params.cutLead + 2 - atom) &&
              j >  (params.leadWidth + 2 - atom)) ||
             (j >= (params.nY + nInFold + params.cutLead - 1 - atom) &&
              j <  (lastY - params.leadWidth - 1 - atom))))
        {
          add = false;
        }

        // cut for X higher lead
        if ((i >= params.nX - 2) &&
            ((j <= (params.nY - params.cutLead + 1 - atom) &&
              j >  (params.leadWidth + 1 - atom)) ||
             (j >= (params.nY + nInFold + params.cutLead - 2 - atom) &&
              j <  (lastY - params.leadWidth - 2 - atom))))
        {
          add = false;
        }

        Vec3S pos = {0., 0., 0.};

        if ((j <= params.nY && atom == 0) ||
            (j <  params.nY && atom == 1))
        {
          pos = genFlatPos(i, j, atom);
        }
        else if ((j <= (params.nY + nInFold) && atom == 0) ||
                 (j <  (params.nY + nInFold) && atom == 1))
        {
          pos = genFold(i, j, positions, R, params);
          // varLog(pos);
        }
        else
        {
          pos = genUpperPos(
            j,
            hasUpperPos ? prevUpperPos : (hasFoldPos ? lastFoldPos : positions.back()),
            R,
            !hasUpperPos);
          prevUpperPos = pos;
          hasUpperPos = true;
        }

        maxZ = std::max(maxZ, pos.z);
        if (add)
        {
          positions.push_back(pos);
        }

        const bool isFoldPosition =
          (atom == 0 && j >= params.nY && j <= params.nY + nInFold) ||
          (atom == 1 && j >= params.nY && j < params.nY + nInFold);
        if (isFoldPosition)
        {
          lastFoldPos = pos;
          hasFoldPos = true;
        }
      }
    }
  }

  return positions;
}

//////////////// local function definitions ////////////////

static Vec3S genFlatPos(const size_t i, const size_t j, const uint8_t atom)
{
  Vec3S pos = {};

  pos.x = Const::atomsArmchair[atom][0];
  pos.y = Const::atomsArmchair[atom][1];

  pos.x += Const::vecsArmchair[0][0] * i;
  pos.y += Const::vecsArmchair[0][1] * i;

  pos.x += Const::vecsArmchair[1][0] * j;
  pos.y += Const::vecsArmchair[1][1] * j;

  pos.x += Const::posOffsetArmchair[0];
  pos.y += Const::posOffsetArmchair[1];

  // correct for rectangle
  pos.x -= 2 * (j / 2) * Const::vecsArmchair[1][0];

  return pos;
}

static Vec3S genUpperPos(
  const size_t j,
  const Vec3S& prevPos,
  const double foldRadius,
  const bool firstInUpper)
{
  const double nextX = prevPos.x +
    Const::vecsArmchair[1][0] * ((j % 2 == 0) ? -1.0 : 1.0);
  const double upperZ = 2.0 * foldRadius;

  if (!firstInUpper)
  {
    return {nextX, prevPos.y - Const::vecsArmchair[1][1], upperZ};
  }

  const double bondLength = std::hypot(
    Const::vecsArmchair[1][0], Const::vecsArmchair[1][1]);
  const double deltaX = nextX - prevPos.x;
  const double deltaZ = upperZ - prevPos.z;
  const double deltaYSquared =
    bondLength * bondLength - deltaX * deltaX - deltaZ * deltaZ;
  const double deltaY = std::sqrt(std::max(0.0, deltaYSquared));

  return
  {
    nextX,
    prevPos.y - deltaY,
    upperZ
  };
}

static Vec3S genFold(
  const size_t i,
  const size_t j,
  const std::vector<Vec3S>& positions,
  const double foldRadius,
  const ParamsS& params)
{
  const Vec3S& prevPos = positions.back();
  const double centerY = genFlatPos(i, params.nY, 0).y;
  const double centerZ = foldRadius;
  const double yzBondLength = Const::vecsArmchair[1][1];

  const double deltaY = prevPos.y - centerY;
  const double deltaZ = prevPos.z - centerZ;
  const double centerDistance = std::hypot(deltaY, deltaZ);
  const double radiusSum = foldRadius + yzBondLength;
  const double radiusDifference = std::abs(foldRadius - yzBondLength);

  if (centerDistance <= radiusDifference || centerDistance >= radiusSum)
  {
    const double radialY = deltaY / centerDistance;
    const double radialZ = deltaZ / centerDistance;
    return
    {
      prevPos.x + Const::vecsArmchair[1][0] * ((j % 2 == 0) ? -1.0 : 1.0),
      centerY + foldRadius * radialY,
      centerZ + foldRadius * radialZ
    };
  }

  const double alongCenters =
    (foldRadius * foldRadius - yzBondLength * yzBondLength +
     centerDistance * centerDistance) / (2.0 * centerDistance);
  const double heightSquared = foldRadius * foldRadius - alongCenters * alongCenters;
  const double height = std::sqrt(std::max(0.0, heightSquared));
  const double baseY = centerY + alongCenters * deltaY / centerDistance;
  const double baseZ = centerZ + alongCenters * deltaZ / centerDistance;

  const double offsetY = -height * deltaZ / centerDistance;
  const double offsetZ = height * deltaY / centerDistance;
  const double firstZ = baseZ + offsetZ;
  // const bool chooseFirstIntersection = (atom == 0);
  const double nextY = baseY + offsetY;// : baseY - offsetY;
  const double nextZ = firstZ; // : secondZ;
  const double nextX = prevPos.x +
    Const::vecsArmchair[1][0] * ((j % 2 == 0) ? -1.0 : 1.0);

  return {nextX, nextY, nextZ};
}


