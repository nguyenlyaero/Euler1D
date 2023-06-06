#include <cmath>
#include <iostream>

#include "idealgaseos.h"

void IdealGasEOS::Initialize(const toml::table& params) {
  this->gamma0 = *params["EOS"]["gamma0"].value<double>();
  this->gamma1 = *params["EOS"]["gamma1"].value<double>();
  this->G0 = 1./(gamma0-1.);
  this->G1 = 1./(gamma1-1.);
  this->M0 = *params["EOS"]["M0"].value<double>();
  this->M1 = *params["EOS"]["M1"].value<double>();

  std::string flux = *params["EOS"]["NumFlux"].value<std::string>();
  if (flux == "CENTRAL") {
    this->numerical_flux = CENTRAL;
  } else if (flux == "PEQ_CENTRAL") {
    this->numerical_flux = PEQ_CENTRAL;
  } else if (flux == "UPWIND_LEFT") {
    this->numerical_flux = UPWIND_LEFT;
  } else if (flux == "PEQ_UPWIND_LEFT") {
    this->numerical_flux = PEQ_UPWIND_LEFT;
  } else {
    std::cout << "Unsupported NumFlux: " << flux << std::endl;
    throw -1;
  }

  this->PrintSelf();
}

void IdealGasEOS::PrintSelf() {
  std::stringstream msg;
  msg << std::endl
      << "------------------------------------" << std::endl
      << "IdealGasEOS::PrintSelf()" << std::endl
      << "------------------------------------" << std::endl;

  msg << "gamma0: " << this->gamma0 << std::endl;
  msg << "gamma1: " << this->gamma1 << std::endl;
  msg << "M0: " << this->M0 << std::endl;
  msg << "M1: " << this->M1 << std::endl;
  if (this->numerical_flux == CENTRAL)
    msg << "Numerical flux: " << "CENTRAL" << std::endl;
  else if (this->numerical_flux == PEQ_CENTRAL)
    msg << "Numerical flux: " << "PEQ_CENTRAL" << std::endl;
  else if (this->numerical_flux == UPWIND_LEFT)
    msg << "Numerical flux: " << "UPWIND_LEFT" << std::endl;
  else if (this->numerical_flux == PEQ_UPWIND_LEFT)
    msg << "Numerical flux: " << "PEQ_UPWIND_LEFT" << std::endl;

  std::cout << msg.str();
}

//-------------------------------------------------------------------------------------

void IdealGasEOS::CalcNumericalFlux(double R0, double R1,
                                    double RU0, double RU1,
                                    double RE0, double RE1,
                                    double RY0, double RY1,
                                    double P0, double P1,
                                    double T0, double T1,
                                    double& Frho0, double& Frhou0,
                                    double& FrhoE0, double& FrhoY0,
                                    double& Frho1, double& Frhou1,
                                    double& FrhoE1, double& FrhoY1,
                                    double BR0, double BR1,
                                    double BRY0, double BRY1,
                                    double& FA_MARK) {
  double Frho, Frhou, FrhoE, FrhoY;

  if (this->numerical_flux == CENTRAL || this->numerical_flux == PEQ_CENTRAL) {
    double G0, G1;
    G0 = (RE0 - 0.5*RU0*RU0/R0)/P0;
    G1 = (RE1 - 0.5*RU1*RU1/R1)/P1;

    double U0 = RU0/R0;
    double U1 = RU1/R1;

    double C0f, Cf, Mf, PIf, Kf, If, Pf;
    if (this->numerical_flux == CENTRAL) {
      C0f = 0.5*(RU0 + RU1);
      Cf = 0.5*(RY0*U0 + RY1*U1);
      Mf = 0.5*(RU0*U0 + RU1*U1);
      PIf = 0.5*(P0 + P1);
      Kf = 0.5*(RU0*0.5*U0*U0 + RU1*0.5*U1*U1);
      If = 0.5*(P0*G0*U0 + P1*G1*U1);
      Pf = 0.5*(U0*P0 + U1*P1);
    } else if (this->numerical_flux == PEQ_CENTRAL) {
      double phi0, phi1;
      this->Y = RY0/R0;
      this->CalcMixingRule();
      double Mbar0 = this->M_bar;

      this->Y = RY1/R1;
      this->CalcMixingRule();
      double Mbar1 = this->M_bar;

      phi0 = Mbar0/Mbar1*R1/R0;
      phi1 = Mbar1/Mbar0*R0/R1;

      double Uf = 0.5*(U0 + U1);
      double rhof = 0.5*(phi0*R0 + phi1*R1);
      double rhoYf = 0.5*(phi0*RY0 + phi1*RY1);

      C0f = rhof * Uf;
      Cf = rhoYf * Uf;
      Mf = rhof * Uf * Uf;
      PIf = 0.5*(P0 + P1);
      Kf = rhof * Uf * 0.5*U0*U1;
      If = 0.5*(P0*G0 + P1*G1) * Uf;
      Pf = 0.5*(U0*P1 + U1*P0);
    }
    Frho = C0f;
    Frhou = Mf + PIf;
    FrhoE = Kf + If + Pf;
    FrhoY = Cf;

    return;
  } else if (this->numerical_flux == UPWIND_LEFT) {
    assert(RU0/R0 > 0. && RU1/R1 > 0.);

    Frho = RU0;
    Frhou = RU0*RU0/R0 + 0.5*(P0+P1);
    FrhoE = RE0*RU0/R0 + 0.5*(P0+P1)*(RU0/R0);
    FrhoY = RY0*RU0/R0;

    return;
  } else if (this->numerical_flux == PEQ_UPWIND_LEFT) {
//    assert(RU0/R0 > 0. && RU1/R1 > 0.);

//    Frho = RU0;
//    Frhou = RU0*RU0/R0 + 0.5*(P0+P1);
//    FrhoE = RE0*RU0/R0 + 0.5*(P0+P1)*(RU0/R0);

    double G0, G1;
    G0 = (RE0 - 0.5*RU0*RU0/R0)/P0;
    G1 = (RE1 - 0.5*RU1*RU1/R1)/P1;

    double U0 = RU0/R0;
    double U1 = RU1/R1;

    double C0f, Cf, Mf, PIf, Kf, If, Pf;
    C0f = 0.5*(RU0 + RU1);
    Cf = 0.5*(RY0*U0 + RY1*U1);
    Mf = 0.5*(RU0*U0 + RU1*U1);
    PIf = 0.5*(P0 + P1);
    Kf = 0.5*(RU0*0.5*U0*U0 + RU1*0.5*U1*U1);
    If = 0.5*(P0*G0*U0 + P1*G1*U1);
    Pf = 0.5*(U0*P0 + U1*P1);

    Frho = C0f;
    Frhou = Mf + PIf;
    FrhoE = Kf + If + Pf;
    FrhoY = Cf;


//    double Re0 = (RE0 - 0.5*RU0*RU0/R0);
//    double Re1 = (RE1 - 0.5*RU1*RU1/R1);
//    double u0 = 0.5*(RU0/R0 + RU1/R1);
//    double numerator = u0 * (Re1 - Re0)
//        - u0 * (BR1*R1 - BR0*R0)
//        - u0 * (BRY1*RY1 - BRY0*RY0)
//        + (BR1 - BR0) * Frho;
//    double denominator = -(BRY1 - BRY0);

//    FrhoY = numerator / denominator;

    double u0 = 0.5*(RU0/R0 + RU1/R1);
    double FR_0 = RU0, FR_1 = RU1;
    double FRU_0 = RU0*RU0/R0 + P0, FRU_1 = RU1*RU1/R1 + P1;
    double FRE_0 = RE0*RU0/R0 + P0*RU0/R0, FRE_1 = RE1*RU1/R1 + P1*RU1/R1;
    double FRY_0 = RY0*RU0/R0, FRY_1 = RY1*RU1/R1;
    double numerator = (FRE_1 - FRE_0) + u0*u0/2.*(FR_1 - FR_0)
        - (BR1*FR_1 - BR0*FR_0)
        + (BR1 - BR0) * Frho
        - u0*(FRU_1 - FRU_0)
        - (BRY1*FRY_1 - BRY0*FRY_0);
    double denominator = -(BRY1 - BRY0);

    FrhoY = numerator / denominator;

    // Limiter
    double r = std::fabs(R0 - R1)/std::max(R0, R1);
    double zero = 0.005;
    double k = 0.005/3;
    FA_MARK = 0.5*(1 + std::tanh((r - zero)/k));

    if (FA_MARK < 1e-2 || std::isnan(FrhoY))
      FrhoY = Cf;
    else
      FrhoY = FA_MARK * FrhoY + (1. - FA_MARK)*Cf;

  } else {
    throw -1;
  }

  Frho0 = Frho; Frho1 = Frho;
  Frhou0 = Frhou; Frhou1 = Frhou;
  FrhoE0 = FrhoE; FrhoE1 = FrhoE;
  FrhoY0 = FrhoY; FrhoY1 = FrhoY;
}

//-------------------------------------------------------------------------------------

void IdealGasEOS::CalcMixingRule() {
  this->M_bar = 1. / (this->Y/this->M0 + (1.-this->Y)/this->M1);
  double cp_bar = this->Y*this->gamma0/(this->gamma0-1.)*this->R0/this->M0
      + (1.-this->Y)*this->gamma1/(this->gamma1-1.)*this->R0/this->M1;
  double cv_bar = this->Y/(this->gamma0-1.)*this->R0/this->M0
      + (1.-this->Y)/(this->gamma1-1.)*this->R0/this->M1;
  this->gamma_bar = cp_bar/cv_bar;

  this->R_bar = this->R0 / this->M_bar;
}

void IdealGasEOS::CalcBetas() {
  this->beta_R = -this->P * this->M_bar*this->M_bar / this->M0/this->M1
      * (this->G0 - this->G1) * this->Y / this->R;
  this->beta_RY = this->P * this->M_bar*this->M_bar / this->M0/this->M1
      * (this->G0 - this->G1) / this->R;
}

void IdealGasEOS::SetMixture_TPY(double T, double P, double Y) {
  this->T = T;
  this->P = P;
  this->Y = Y;

  this->CalcMixingRule();

  this->R = this->P / this->T / this->R_bar;
  this->E = this->P / this->R / (this->gamma_bar - 1);

  this->SoS = std::sqrt(this->gamma_bar * this->P / this->R);

  this->CalcBetas();
}

void IdealGasEOS::SetMixture_TRY(double T, double R, double Y) {
  this->T = T;
  this->R = R;
  this->Y = Y;

  this->CalcMixingRule();

  this->P = this->R * this->R_bar * this->T;
  this->E = this->P / this->R / (this->gamma_bar - 1);

  this->SoS = std::sqrt(this->gamma_bar * this->P / this->R);

  this->CalcBetas();
}

void IdealGasEOS::SetMixture_ERY(double E, double R, double Y, double Tguess) {
  this->E = E;
  this->R = R;
  this->Y = Y;

  this->CalcMixingRule();

  this->P = this->E * this->R * (this->gamma_bar - 1);
  this->T = this->P / this->R / this->R_bar;

  this->SoS = std::sqrt(this->gamma_bar * this->P / this->R);

  this->CalcBetas();
}

void IdealGasEOS::SetMixture_PRY(double P, double R, double Y) {
  this->P = P;
  this->R = R;
  this->Y = Y;

  this->CalcMixingRule();

  this->E = this->P / this->R / (this->gamma_bar - 1);
  this->T = this->P / this->R / this->R_bar;

  this->SoS = std::sqrt(this->gamma_bar * this->P / this->R);

  this->CalcBetas();
}


