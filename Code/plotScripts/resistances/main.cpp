#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fstream>
#include <cassert>

#include "resistances.h"

void handleArgs(int argc, char* argv[], std::string& directory, LeadInfoS& leadInfo);

int main(int argc, char* argv[])
{
  std::string directory = "./results/";

  LeadInfoS leadInfo {
    .currentLeadFrom = 3,
    .currentLeadTo = 4,
    .voltageLeadHigh = 2,
    .voltageLeadLow = 1
  };

  handleArgs(argc, argv, directory, leadInfo);

  std::cout << "Reading results form directory " << directory << '\n';
  const std::filesystem::path dirs(directory + "dirs/");

  size_t numOfDirs = getNumOfSubDirs(dirs);
  std::vector<double> results;
  results.reserve(numOfDirs);

  // remapping is done to keep way of calculating of resistanced based on Datta (eradicate row and column 4)
  const bool needsRemap = leadDesired != leadInfo;
  std::cout << "needsRemap: " << needsRemap << "\n";

  for (auto& singleResultDir: std::filesystem::directory_iterator(dirs))
  {
    std::filesystem::path filePath(singleResultDir.path()/"Transmissions.csv");
    dmsg("Reading results form file " << filePath);

    double T[16] = {};

    // T_to_from
    double R = 0;
    if (!readMatrix(filePath, T, leadInfo, needsRemap))
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
};

void handleArgs(int argc, char* argv[], std::string& directory, LeadInfoS& leadInfo)
{
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

  bool valid = true;
  if (argc == 6)
  {
    leadInfo.currentLeadFrom = std::atoi(argv[2]);
    leadInfo.currentLeadTo = std::atoi(argv[3]);
    leadInfo.voltageLeadHigh = std::atoi(argv[4]);
    leadInfo.voltageLeadLow = std::atoi(argv[5]);
    valid &= (leadInfo.currentLeadFrom > 0 && leadInfo.currentLeadFrom <= 4);
    valid &= (leadInfo.currentLeadTo > 0 && leadInfo.currentLeadTo <= 4);
    valid &= (leadInfo.voltageLeadHigh > 0 && leadInfo.voltageLeadHigh <= 4);
    valid &= (leadInfo.voltageLeadLow > 0 && leadInfo.voltageLeadLow <= 4);
    valid &= (leadInfo.currentLeadFrom != leadInfo.currentLeadTo);
    valid &= (leadInfo.currentLeadFrom != leadInfo.voltageLeadLow);
    valid &= (leadInfo.currentLeadFrom != leadInfo.voltageLeadHigh);
    valid &= (leadInfo.currentLeadTo   != leadInfo.voltageLeadLow);
    valid &= (leadInfo.currentLeadTo   != leadInfo.voltageLeadHigh);
    valid &= (leadInfo.voltageLeadHigh != leadInfo.voltageLeadLow);
  }
  if (!valid)
  {
    std::cout << "args not valid\n";
    std::exit(1);
  }
}

