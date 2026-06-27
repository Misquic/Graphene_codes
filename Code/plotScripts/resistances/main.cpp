#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fstream>
#include <cassert>

#include "resistances.h"

int main(int argc, char* argv[])
{
  std::string directory = "./results/";

  if (argc >= 2)
  {
    directory = argv[1];

    if (directory == "help" || directory == "-h" || !(argc == 2 || argc == 6))
    {
      std::printf("directory currentLeadFrom currentLeadTo voltageLeadHigh voltageLeadLow\n");
      std::exit(0);
    }

    if (directory.back() != '/')
      {
        directory += '/';
      }
  }

  leadInfoS leadInfo {
    .currentLeadFrom = 1,
    .currentLeadTo = 2,
    .voltageLeadHigh = 3,
    .voltageLeadLow = 4
  };

  if (argc == 6)
  {
    leadInfo.currentLeadFrom = std::atoi(argv[2]);
    leadInfo.currentLeadTo = std::atoi(argv[3]);
    leadInfo.voltageLeadHigh = std::atoi(argv[4]);
    leadInfo.voltageLeadLow = std::atoi(argv[5]);
    assert(leadInfo.currentLeadFrom > 0 && leadInfo.currentLeadFrom <= 4);
    assert(leadInfo.currentLeadTo > 0 && leadInfo.currentLeadTo <= 4);
    assert(leadInfo.voltageLeadHigh > 0 && leadInfo.voltageLeadHigh <= 4);
    assert(leadInfo.voltageLeadLow > 0 && leadInfo.voltageLeadLow <= 4);
  }

  std::cout << "Reading results form directory " << directory << '\n';
  const std::filesystem::path dirs(directory + "dirs/");

  size_t numOfDirs = getNumOfSubDirs(dirs);
  std::vector<double> results;
  results.reserve(numOfDirs);

  for (auto& singleResultDir: std::filesystem::directory_iterator(dirs))
  {
    std::filesystem::path filePath(singleResultDir.path()/"Transmissions.csv");
    dmsg("Reading results form file " << filePath);

    double T[16] = {};

    double R = 0;
    if (!readMatrix(filePath, T, leadInfo))
    {
      std::cerr << "Matrix not read properly, setting R = -1\n";
      R = -1;
    }
    else
    {
      R = calculateR(T);
    }

    results.push_back(R);
    dmsg(R);
  }

  bool valid = true;
  std::ofstream file = openOFile(directory + "R.dat", valid);
  if (!valid) exit(1);

  int i = 0;
  for (auto& singleResultDir: std::filesystem::directory_iterator(dirs))
  {
    file << singleResultDir.path().filename() << ',' << results[i] << '\n';
    i++;
  }
  // for (size_t i = 0; i < results.size(); i++)
  // {
    // file << results[i] << '\n';
  // }

};
