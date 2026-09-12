#pragma once

#include "Utilities.h"
#include "Constants.h"

typedef struct ParamsS
{
  size_t nX = 10;
  size_t nY = 10;
  double foldRadius = 10.f;
  size_t cutLead = 5;
  size_t leadWidth = 5;
  std::string saveFileName = "./results/system.csv";
} ParamsS;

typedef struct Vec3S
{
  double x;
  double y;
  double z;
} Vec3S;

inline std::ostream& operator<<(std::ostream& out, const Vec3S vec);

void readArgs(const int argc, const char* const* const argv, ParamsS& params);

std::vector<Vec3S>& generatePositions(const ParamsS& params);

void printPositions(const std::vector<Vec3S>& positions);

void savePositions(const std::string_view name, const std::vector<Vec3S>& positions);

void makeParamsRight(ParamsS& params);

//////////////////// inline functions ///////////////////

inline Vec3S operator+(const Vec3S vec1, const Vec3S vec2)
{
  return {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z};
}

inline Vec3S operator-(const Vec3S vec1, const Vec3S vec2)
{
  return {vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z};
}

inline void operator+=(Vec3S& vec1, const Vec3S vec2)
{
  vec1.x += vec2.x;
  vec1.y += vec2.y;
  vec1.z += vec2.z;
}

inline void operator-=(Vec3S& vec1, const Vec3S vec2)
{
  vec1.x -= vec2.x;
  vec1.y -= vec2.y;
  vec1.z -= vec2.z;
}

inline double operator*(const Vec3S vec1, const Vec3S vec2)
{
  return {vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z};
}

inline double vecLen(const Vec3S vec1)
{
  double len = std::sqrt(vec1 * vec1);
  dmsg("vec = " << vec1 << " len = " << len);
  return len;
}

inline constexpr double deg2rad(double deg)
{
  return deg * M_PI / 180.f;
}

inline constexpr double rad2deg(double rad)
{
  return rad * 180.f / M_PI;
}

inline double getAngleBetweenVecsInRad(const Vec3S vec1, const Vec3S vec2)
{
  const double dot = vec1 * vec2;
  const double len1len2 = vecLen(vec1) * vecLen(vec2);

  return std::acos(std::clamp(dot / len1len2, -1.0, 1.0));
}

inline Vec3S rotYZplane(const Vec3S vec, double radians)
{
  float sin = std::sin(radians);
  float cos = std::cos(radians);

  return
  {
    vec.x,
    vec.y * cos - vec.z * sin,
    vec.y * sin + vec.z * cos,
  };
}

inline std::ostream& operator<<(std::ostream& out, const Vec3S vec)
{
  out << vec.x << ", " << vec.y << ", " << vec.z;
  return out;
}
