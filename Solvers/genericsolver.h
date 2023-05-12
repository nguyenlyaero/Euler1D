#ifndef GENERICSOLVER_H
#define GENERICSOLVER_H

#include <vector>
#include "toml++/toml.hpp"

#include "Interface.h"

#include <GridOps/grid1d.h>
#include <EOS/genericeos.h>

enum DtModes { CONSTANT_DT, CONSTANT_CFL };
enum StoppingMode { NSTEPS, FINAL_TIME};

class GenericSolver {
public:
  GenericSolver() = default;
  virtual ~GenericSolver();

  virtual void Initialize(const toml::table& params); // TODO

  virtual void PrintSelf();

  virtual void RunSolver();

protected:
  virtual void InitializeCase(const toml::table& params);

  void CalcDt();
  void PrintTimeStepInfo();
  bool DoneSolver();

  void SolveRKSubStep(int rk_step);
  void CalcRKRhs(int rk_step, double rk_time);

  void UpdateConservatives(int rk_step);
  void UpdatePrimitives(int rk_step);
  void UpdateFluxes(int rk_step);

  void Output();
  static void DumpScalarRange(double* scalar, int n, std::string name);


  // RK coefficients
  constexpr static const int N_RK_STEP = 3;
  constexpr static const double RK_COEFF[3][3] =
    {{1.0, 0., 0.}, {-3. / 4., 1. / 4., 0.}, {-1. / 12., -1. / 12., 8. / 12.}};
  constexpr static const double RK_TIME_COEFF[3] = {1.0, 0.5, 0.};

  Grid1D grid;
  int ncv, nfa;
  const std::pair<int,int>* cvofa;
  const double* xcv;
  const double* cv_volume;

  GenericEOS* mixture;

  constexpr static const int num_conservatives = 4;
  double* R;
  double* RU;
  double* RE;
  double* RY;

  double* P;
  double* T;
  double* SoS;
  double* CFL;

  double* beta_R;
  double* beta_RY;

  // Fluxes
  double* FR;
  double* FRU;
  double* FRE;
  double* FRY;

  double* FA_MARK;

  double* FR_CV;
  double* FRU_CV;
  double* FRE_CV;
  double* FRY_CV;

  // Arrays
  double** conservatives;
  double*** rhs_array;

  // Time stepping
  int step;
  double time;

  DtModes dt_mode;
  double CFL_input;
  double dt_input;
  double dt_actual;

  // Check done
  StoppingMode stop;
  double final_time;
  int total_steps;

  int write_interval;
  std::vector<std::string> output_vars;
  std::vector<std::string> output_vars_fa;
};

#endif // GENERICSOLVER_H
