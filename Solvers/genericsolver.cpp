#include "genericsolver.h"

#include <algorithm>
#include <iomanip>
#include <fstream>

GenericSolver::~GenericSolver() {
  delete[] this->conservatives;

  for (int i = 0; i < this->N_RK_STEP; ++i)
    delete[] this->rhs_array[i];
  delete[] this->rhs_array;

  delete mixture;
}

//-------------------------------------------------------------------------------------

void GenericSolver::Initialize(const toml::table& params) {
  this->grid.Initialize(params);
  this->mixture = MakeEOS(params);

  std::string timestep = *params["Solve"]["TimeStepping"].value<std::string>();
  if (timestep == "CONSTANT_DT") {
    this->dt_mode = CONSTANT_DT;
    this->dt_input = *params["Solve"]["dt"].value<double>();
  } else if (timestep == "CONSTANT_CFL") {
    this->dt_mode = CONSTANT_CFL;
    this->CFL_input = *params["Solve"]["CFL"].value<double>();
  } else {
    std::cout << "Unsupported TimeStepping: " << timestep << std::endl;
    throw -1;
  }

  std::string stopping = *params["Solve"]["Stopping"].value<std::string>();
  if (stopping == "NSTEPS") {
    this->stop = NSTEPS;
    this->total_steps = *params["Solve"]["num_steps"].value<int>();
  } else if (stopping == "FINAL_TIME") {
    this->stop = FINAL_TIME;
    this->final_time = *params["Solve"]["final_time"].value<double>();
  } else {
    std::cout << "Unsupported Stopping: " << stopping << std::endl;
    throw -1;
  }

  {
    const auto& stringList = *params["Solve"]["Output"].as_array();
    for (const auto& element : stringList) {
      this->output_vars.push_back(*element.value<std::string>());
    }
  }
  {
    const auto& stringList = *params["Solve"]["Output_FA"].as_array();
    for (const auto& element : stringList) {
      this->output_vars_fa.push_back(*element.value<std::string>());
    }
  }
  this->write_interval = *params["Solve"]["write_interval"].value<int>();

  this->PrintSelf();

  this->ncv = this->grid.GetNCV();
  this->nfa = this->grid.GetNFA();
  this->cvofa = this->grid.GetCVOFA();
  this->xcv = this->grid.GetXCV();
  this->cv_volume = this->grid.GetCVVolume();

  // Register variables
  this->R = this->grid.RegisterVariable("RHO", CV_DATA);
  this->RU = this->grid.RegisterVariable("RHOU", CV_DATA);
  this->RE = this->grid.RegisterVariable("RHOE", CV_DATA);
  this->RY = this->grid.RegisterVariable("RHOY", CV_DATA);

  this->P = this->grid.RegisterVariable("P", CV_DATA);
  this->T = this->grid.RegisterVariable("T", CV_DATA);
  this->SoS = this->grid.RegisterVariable("SOS", CV_DATA);
  this->CFL = this->grid.RegisterVariable("CFL", CV_DATA);

  this->beta_R = this->grid.RegisterVariable("beta_R", CV_DATA);
  this->beta_RY = this->grid.RegisterVariable("beta_RY", CV_DATA);

  this->conservatives = new double*[this->num_conservatives];
  this->conservatives[0] = this->R;
  this->conservatives[1] = this->RU;
  this->conservatives[2] = this->RE;
  this->conservatives[3] = this->RY;

  this->rhs_array = new double**[this->N_RK_STEP];
  for (int rk_step = 0; rk_step < this->N_RK_STEP; ++rk_step) {
    this->rhs_array[rk_step] = new double*[this->num_conservatives];
    this->rhs_array[rk_step][0] = this->grid.RegisterVariable("RHS_RHO_"+std::to_string(rk_step), CV_DATA);
    this->rhs_array[rk_step][1] = this->grid.RegisterVariable("RHS_RHOU_"+std::to_string(rk_step), CV_DATA);
    this->rhs_array[rk_step][2] = this->grid.RegisterVariable("RHS_RHOE_"+std::to_string(rk_step), CV_DATA);
    this->rhs_array[rk_step][3] = this->grid.RegisterVariable("RHS_RHOY_"+std::to_string(rk_step), CV_DATA);
  }

  this->FR0 = this->grid.RegisterVariable("FR0", FA_DATA);
  this->FRU0 = this->grid.RegisterVariable("FRU0", FA_DATA);
  this->FRE0 = this->grid.RegisterVariable("FRE0", FA_DATA);
  this->FRY0 = this->grid.RegisterVariable("FRY0", FA_DATA);

  this->FR1 = this->grid.RegisterVariable("FR1", FA_DATA);
  this->FRU1 = this->grid.RegisterVariable("FRU1", FA_DATA);
  this->FRE1 = this->grid.RegisterVariable("FRE1", FA_DATA);
  this->FRY1 = this->grid.RegisterVariable("FRY1", FA_DATA);

  this->FA_MARK = this->grid.RegisterVariable("FA_MARK", FA_DATA);

  this->FR_CV = this->grid.RegisterVariable("FR_CV", CV_DATA);
  this->FRU_CV = this->grid.RegisterVariable("FRU_CV", CV_DATA);
  this->FRE_CV = this->grid.RegisterVariable("FRE_CV", CV_DATA);
  this->FRY_CV = this->grid.RegisterVariable("FRY_CV", CV_DATA);


  this->PrintSelf();

  // TODO: Initialize case
  this->step = 0;
  this->time = 0.;

  this->InitializeCase(params);
}

void GenericSolver::InitializeCase(const toml::table& params) {
  std::string init = *params["IC"]["IC_type"].value<std::string>();
  if (init == "Kawai_Ramp") {
    double u0 = *params["IC"]["u0"].value<double>();
    double P0 = *params["IC"]["P0"].value<double>();
    double xc = *params["IC"]["xc"].value<double>();
    double rc = *params["IC"]["rc"].value<double>();
    double k = *params["IC"]["k"].value<double>();
    double wY0 = *params["IC"]["wY0"].value<double>();
    double wY1 = *params["IC"]["wY1"].value<double>();

    std::cout << "------------------------------------" << std::endl;
    std::cout << "Initialize Case: " << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "u0: " << u0 << ", P0: " << P0 << std::endl;
    std::cout << "xc: " << xc << ", rc: " << rc << ", k: " << k << std::endl;
    std::cout << "wY0: " << wY0 << ", wY1: " << wY1 << std::endl;

    for (int icv = 0; icv < this->ncv; ++icv) {
      double r = std::fabs(this->xcv[icv] - xc);
      double RY0 = wY0/2.*(1. - std::tanh(k*(r - rc)));
      double RY1 = wY1/2.*(1. + std::tanh(k*(r - rc)));

      this->R[icv] = RY0 + RY1;
      this->RY[icv] = RY0;
      this->RU[icv] = this->R[icv] * u0;
      this->mixture->SetMixture_PRY(P0, this->R[icv], this->RY[icv]/this->R[icv]);
      this->RE[icv] = this->R[icv]*this->mixture->GetE()
          + 0.5*this->R[icv]*u0*u0;
    }
  } else if (init == "TEST_TP_Linear") {
    double T0 = *params["IC"]["T0"].value<double>();
    double T1 = *params["IC"]["T1"].value<double>();
    assert(T0 < T1);
    double x0 = *params["IC"]["x0"].value<double>();
    double x1 = *params["IC"]["x1"].value<double>();
    assert(x0 < x1);
    double P0 = *params["IC"]["P0"].value<double>();

    std::cout << "------------------------------------" << std::endl;
    std::cout << "Initialize Case: " << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "u0: " << 0 << ", P: " << P0 << std::endl;
    std::cout << "T0: " << T0 << ", T1: " << T1 << std::endl;
    std::cout << "x0: " << x0 << ", x1: " << x1 << std::endl;

    for (int icv = 0; icv < this->ncv; ++icv) {
      double Tin = std::min(std::max(T0, T0 + (T1 - T0) / (x1 - x0) * (xcv[icv] - x0)), T1);

      this->mixture->SetMixture_TPY(Tin, P0, 0.);

      this->R[icv] = this->mixture->GetR();
      this->RY[icv] = 0.;
      this->RU[icv] = 0.;
      this->RE[icv] = this->R[icv]*this->mixture->GetE();
    }
  } else if (init == "TP_Ramp") {
    double T0 = *params["IC"]["T0"].value<double>();
    double T1 = *params["IC"]["T1"].value<double>();
    double P0 = *params["IC"]["P0"].value<double>();
    double u0 = *params["IC"]["u0"].value<double>();
    double xc = *params["IC"]["xc"].value<double>();
    double rc = *params["IC"]["rc"].value<double>();
    double k = *params["IC"]["k"].value<double>();

    std::cout << "------------------------------------" << std::endl;
    std::cout << "Initialize Case: " << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "u0: " << u0 << ", P: " << P0 << std::endl;
    std::cout << "T0: " << T0 << ", T1: " << T1 << std::endl;
    std::cout << "xc: " << xc << ", r1: " << rc << ", k: " << k <<  std::endl;

    for (int icv = 0; icv < this->ncv; ++icv) {
      double r = std::fabs(this->xcv[icv] - xc);
      double Tin = T0 + 0.5 * (T1 - T0) * (1. + std::tanh(k*(r - rc)));
      this->mixture->SetMixture_TPY(Tin, P0, 0.);

      this->R[icv] = this->mixture->GetR();
      this->RY[icv] = 0.;
      this->RU[icv] = this->R[icv] * u0;
      this->RE[icv] = this->R[icv]*this->mixture->GetE() + 0.5*this->R[icv]*u0*u0;;
    }
  } else {
    std::cout << "Unsupported Initialize Case: " << init << std::endl;
    throw -1;
  }
  this->UpdatePrimitives(0);
  this->UpdateFluxes(0);
}

void GenericSolver::PrintSelf() {
  std::stringstream msg;
  msg << std::endl
      << "------------------------------------" << std::endl
      << "GenericSolver::PrintSelf()" << std::endl
      << "------------------------------------" << std::endl;

  if (this->dt_mode == CONSTANT_DT) {
    msg << "TimeStepping: " << "CONSTANT_DT" << std::endl;
    msg << "dt: " << this->dt_input << std::endl;
  } else if (this->dt_mode == CONSTANT_CFL) {
    msg << "TimeStepping: " << "CONSTANT_CFL" << std::endl;
    msg << "CFL: " << this->CFL_input << std::endl;
  }

  if (this->stop == FINAL_TIME) {
    msg << "Stopping: " << "FINAL_TIME" << std::endl;
    msg << "final_time: " << this->final_time << std::endl;
  } else if (this->stop == NSTEPS) {
    msg << "Stopping: " << "NSTEPS" << std::endl;
    msg << "num_steps: " << this->total_steps << std::endl;
  }

  std::cout << "Output: ";
  for (auto str : this->output_vars) {
    std::cout << str << " ";
  }
  std::cout << std::endl;
  std::cout << "write_interval: " << this->write_interval << std::endl;

  std::cout << msg.str();
}

//-------------------------------------------------------------------------------------

void GenericSolver::CalcDt() {
  if (this->dt_mode == CONSTANT_DT) {
    this->dt_actual = dt_input;
  } else if (this->dt_mode == CONSTANT_CFL) {
    double max_sr = 0.;
    for (int ifa = 1; ifa < this->nfa; ++ifa) {
      int icv0 = this->cvofa[ifa].first;
      int icv1 = this->cvofa[ifa].second;
      double h = this->xcv[icv1] - this->xcv[icv0];
      double sos = 0.5*(this->SoS[icv0] + this->SoS[icv1]);
      max_sr = std::max(max_sr, sos/h);
    }
    this->dt_actual = CFL_input / max_sr;
  } else {
    std::cout << "Unknown dt policy" << std::endl;
    throw -1;
  }
}

void GenericSolver::PrintTimeStepInfo() {
  std::cout << "\n----------------------------------------------------------"
            << std::endl;
  std::cout << " starting step: " << this->step
            << " time: " << this->time << " dt: " << this->dt_actual
            << std::endl;
  std::cout << "----------------------------------------------------------"
            << std::endl;
}

bool GenericSolver::DoneSolver() {
  if (this->stop == NSTEPS)
    return (this->step >= this->total_steps);
  else if (this->stop == FINAL_TIME)
    return (this->time >= this->final_time);
  else {
    std::cout << "Unknown stopping criterion" << std::endl;
    throw -1;
  }
}

//-------------------------------------------------------------------------------------

void GenericSolver::CalcRKRhs(int rk_step, double rk_time) {
  // RHS Array pointers
  double* rho_rhs = this->rhs_array[rk_step][0];
  double* rhou_rhs = this->rhs_array[rk_step][1];
  double* rhoE_rhs = this->rhs_array[rk_step][2];
  double* rhoY_rhs = this->rhs_array[rk_step][3];

  // Make them zero
  for (int icv = 0; icv < this->ncv; ++icv) {
    rho_rhs[icv] = 0.;
    rhou_rhs[icv] = 0.;
    rhoE_rhs[icv] = 0.;
    rhoY_rhs[icv] = 0.;
  }

  // Loop over internal faces
  for (int ifa = 0; ifa < this->nfa; ++ifa) {

    int icv0 = this->cvofa[ifa].first;
    int icv1 = this->cvofa[ifa].second;

    rho_rhs[icv0] -= FR0[ifa];
    rhou_rhs[icv0] -= FRU0[ifa];
    rhoE_rhs[icv0] -= FRE0[ifa];
    rhoY_rhs[icv0] -= FRY0[ifa];

    rho_rhs[icv1] += FR1[ifa];
    rhou_rhs[icv1] += FRU1[ifa];
    rhoE_rhs[icv1] += FRE1[ifa];
    rhoY_rhs[icv1] += FRY1[ifa];
  }
}

void GenericSolver::SolveRKSubStep(int rk_step) {
  double rk_time = this->time - RK_TIME_COEFF[rk_step]*this->dt_actual;

  this->CalcRKRhs(rk_step, rk_time);

  // Update conservatives
  this->UpdateConservatives(rk_step);

  // Update primitives
  this->UpdatePrimitives(rk_step);

  // Update fluxes
  this->UpdateFluxes(rk_step);
}

void GenericSolver::UpdateConservatives(int rk_step) {
  for (int isca = 0; isca < this->num_conservatives; ++isca) {
    for (int icv = 0; icv < this->ncv; ++icv) {
      double dt_over_vol = this->dt_actual / this->cv_volume[icv];
      this->rhs_array[rk_step][isca][icv] *= dt_over_vol;

      for (int rk_coeff_index = 0; rk_coeff_index < N_RK_STEP; ++rk_coeff_index) {
        this->conservatives[isca][icv] += RK_COEFF[rk_step][rk_coeff_index] * this->rhs_array[rk_coeff_index][isca][icv];
      }
    }
  }
}

void GenericSolver::UpdatePrimitives(int rk_step) {
  for (int icv = 0; icv < this->ncv; ++icv) {
    double e = RE[icv]/R[icv] - 0.5*RU[icv]*RU[icv]/R[icv]/R[icv];
    this->mixture->SetMixture_ERY(e, R[icv], RY[icv]/R[icv]);

    this->P[icv] = mixture->GetP();
    this->T[icv] = mixture->GetT();
    this->SoS[icv] = mixture->GetSoS();
    this->beta_R[icv] = mixture->GetBeta_R();
    this->beta_RY[icv] = mixture->GetBeta_RY();
  }
}

void GenericSolver::UpdateFluxes(int rk_step) {
  // Loop over internal faces
  for (int ifa = 0; ifa < this->nfa; ++ifa) {
    int icv0 = this->cvofa[ifa].first;
    int icv1 = this->cvofa[ifa].second;

    double rho0 = this->R[icv0];
    double rho1 = this->R[icv1];

    double rhou0 = this->RU[icv0];
    double rhou1 = this->RU[icv1];

    double rhoE0 = this->RE[icv0];
    double rhoE1 = this->RE[icv1];

    double rhoY0 = this->RY[icv0];
    double rhoY1 = this->RY[icv1];

    double P0 = this->P[icv0];
    double P1 = this->P[icv1];

    double T0 = this->T[icv0];
    double T1 = this->T[icv1];

    double BR0 = this->beta_R[icv0];
    double BR1 = this->beta_R[icv1];
    double BRY0 = this->beta_RY[icv0];
    double BRY1 = this->beta_RY[icv1];

    double Frho0 = 0., Frhou0 = 0., FrhoE0 = 0., FrhoY0 = 0.;
    double Frho1 = 0., Frhou1 = 0., FrhoE1 = 0., FrhoY1 = 0.;

    this->mixture->CalcNumericalFlux(rho0, rho1,
                                     rhou0, rhou1,
                                     rhoE0, rhoE1,
                                     rhoY0, rhoY1,
                                     P0, P1,
                                     T0, T1,
                                     Frho0, Frhou0,
                                     FrhoE0, FrhoY0,
                                     Frho1, Frhou1,
                                     FrhoE1, FrhoY1,
                                     BR0, BR1,
                                     BRY0, BRY1, FA_MARK[ifa]);
    this->FR0[ifa] = Frho0;
    this->FRU0[ifa] = Frhou0;
    this->FRE0[ifa] = FrhoE0;
    this->FRY0[ifa] = FrhoY0;

    this->FR1[ifa] = Frho1;
    this->FRU1[ifa] = Frhou1;
    this->FRE1[ifa] = FrhoE1;
    this->FRY1[ifa] = FrhoY1;

    for (int icv : {icv0,icv1}) {
      FR_CV[icv] = RU[icv];
      FRU_CV[icv] = RU[icv]*RU[icv]/R[icv] + P[icv];
      FRE_CV[icv] = RE[icv]*RU[icv]/R[icv] + P[icv]*RU[icv]/R[icv];
      FRY_CV[icv] = RY[icv]*RU[icv]/R[icv];
    }
  }
}



void GenericSolver::RunSolver() {
  bool done = this->DoneSolver();

  this->Output();

  while (!done) {
    ++this->step;
    this->CalcDt();
    this->time += this->dt_actual;
    this->PrintTimeStepInfo();

    for (int i = 0; i < N_RK_STEP; ++i) {
      this->SolveRKSubStep(i);
    }

    this->Output();
    done = this->DoneSolver();
  }
}

//-------------------------------------------------------------------------------------

void GenericSolver::DumpScalarRange(double* scalar, int n, std::string name) {
  double min = *std::min_element(scalar, scalar+n);
  double max = *std::max_element(scalar, scalar+n);

  std::cout << " > dumpScalarRange: " << name << ", " << min << ":" << max << std::endl;
}

void GenericSolver::Output() {
  // Dump scalar range
  GenericSolver::DumpScalarRange(this->R, this->ncv, "R");
  GenericSolver::DumpScalarRange(this->RU, this->ncv, "RU");
  GenericSolver::DumpScalarRange(this->RE, this->ncv, "RE");
  GenericSolver::DumpScalarRange(this->RY, this->ncv, "RY");

  GenericSolver::DumpScalarRange(this->P, this->ncv, "P");
  GenericSolver::DumpScalarRange(this->T, this->ncv, "T");
  GenericSolver::DumpScalarRange(this->SoS, this->ncv, "SoS");

  // Write output
  {
    if (this->step % this->write_interval == 0) {
      std::ofstream fs;
      fs.open("Line."+std::to_string(this->step));
      // Preamble
      fs << "t";
      for (auto str : this->output_vars) fs << "," << str;
      fs << std::endl;

      for (int icv = 0; icv < this->ncv; ++icv) {
        fs << std::setprecision(9) << this->time;
        for(auto str : this->output_vars)
          fs << "," << this->grid.GetCVVar(str)[icv] ;
        fs << std::endl;
      }

      fs.close();
    }
  }
  {
    if (this->step % this->write_interval == 0) {
      std::ofstream fs;
      fs.open("LineFA."+std::to_string(this->step));
      // Preamble
      fs << "t";
      for (auto str : this->output_vars_fa) fs << "," << str;
      fs << std::endl;

      for (int ifa = 0; ifa < this->nfa; ++ifa) {
        fs << std::setprecision(9) << this->time;
        for(auto str : this->output_vars_fa)
          fs << "," << this->grid.GetFAVar(str)[ifa] ;
        fs << std::endl;
      }

      fs.close();
    }
  }
}




