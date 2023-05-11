#ifndef GRID1D_H
#define GRID1D_H

#include <vector>
#include <Eigen/Dense>
#include <unordered_map>
#include <iostream>

#include "toml++/toml.hpp"

#define VALUE_DATA 0
#define CV_DATA    1
#define FA_DATA    2

class Grid1D {
public:
  Grid1D() = default;

  void Initialize(const toml::table& params);

  double* RegisterVariable(std::string name, int data_type);

  double* GetCVVar(std::string name);
  double* GetFAVar(std::string name);

  int GetNCV() const { return this->ncv; }
  int GetNFA() const { return this->nfa; }
  const std::pair<int,int>* GetCVOFA() const { return this->cvofa.data(); }
  const double* GetXCV() const { return this->xcv.data(); }
  const double* GetCVVolume() const { return this->cv_volume.data(); }

  void PrintSelf();
private:
  int ncv, nfa, nno;
  std::vector<double> xcv;
  std::vector<double> xno;
  std::vector<double> xfa;
  std::vector<std::pair<int,int>> cvofa;
  std::vector<double> cv_volume;

  std::vector<std::vector<double>> cv_variables;
  std::unordered_map<std::string, int> cv_variables_indices;
  std::vector<std::vector<double>> fa_variables;
  std::unordered_map<std::string, int> fa_variables_indices;
};

#endif // GRID1D_H
