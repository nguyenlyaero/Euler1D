#include "grid1d.h"

void Grid1D::Initialize(const toml::table& params) {
  int ncv = *params["Grid"]["ncv"].value<int>();
  double xmin = *params["Grid"]["xmin"].value<double>();
  double xmax = *params["Grid"]["xmax"].value<double>();

  this->ncv = ncv;
  this->nfa = ncv;
  this->nno = ncv+1;
  Eigen::VectorXd xno = Eigen::VectorXd::LinSpaced(ncv+1, xmin, xmax);

  this->xno.resize(this->nno);
  for (int ino = 0; ino < this->nno; ++ino)
    this->xno[ino] = xno[ino];

  this->xcv.resize(this->ncv);
  this->cv_volume.resize(this->ncv);
  for (int icv = 0; icv < this->ncv; ++icv) {
    this->xcv[icv] = 0.5*(xno[icv]+xno[icv+1]);
    this->cv_volume[icv] = this->xno[icv+1] - this->xno[icv];
  }

  this->xfa.resize(this->ncv);
  for (int ifa = 0; ifa < this->nfa; ++ifa)
    this->xfa[ifa] = xno[ifa];

  this->cvofa.resize(this->nfa);
  this->cvofa[0] = {this->ncv-1, 0};
  for (int ifa = 1; ifa < this->nfa; ++ifa)
    this->cvofa[ifa] = {ifa-1, ifa};

  this->PrintSelf();
}

void Grid1D::PrintSelf() {
  std::stringstream msg;
  msg << std::endl
      << "------------------------------------" << std::endl
      << "Grid1D::PrintSelf()" << std::endl
      << "------------------------------------" << std::endl;

  msg << "ncv: " << this->ncv << std::endl;
  msg << "nfa: " << this->nfa << std::endl;
  msg << "nno: " << this->nno << std::endl;


  std::cout << msg.str();
}

//-------------------------------------------------------------------------------------

double* Grid1D::RegisterVariable(std::string name, int data_type) {
  std::unordered_map<std::string, int>* map;
  std::vector<std::vector<double>>* vector;
  int n;
  if (data_type == CV_DATA) {
    map = &this->cv_variables_indices;
    vector = &this->cv_variables;
    n = this->ncv;
  } else if (data_type == FA_DATA) {
    map = &this->fa_variables_indices;
    vector = &this->fa_variables;
    n = this->nfa;
  } else {
    std::cout << "Incorrect data type" << std::endl;
    throw -1;
  }

  if (map->find(name) != map->end()) {
    std::cout << "Data name: " << name << " already assigned" << std::endl;
    throw -1;
  }

  map->emplace(name, vector->size());
  vector->push_back(std::vector<double> (n, 0.));

  return (*vector)[vector->size()-1].data();
}

//-------------------------------------------------------------------------------------

double* Grid1D::GetCVVar(std::string name) {
  if (this->cv_variables_indices.find(name) == this->cv_variables_indices.end()) {
    std::cout << "Name " << name << " not found" << std::endl;
    throw -1;
  }
  return this->cv_variables[this->cv_variables_indices[name]].data();
}

double* Grid1D::GetFAVar(std::string name) {
  if (this->fa_variables_indices.find(name) == this->fa_variables_indices.end()) {
    std::cout << "Name " << name << " not found" << std::endl;
    throw -1;
  }
  return this->fa_variables[this->fa_variables_indices[name]].data();
}
