#include <iostream>

#include "vanderwaalseos.h"
#include "utils.h"

void VanDerWaalsEOS::Initialize(const toml::table &params) {
  this->M = *params["EOS"]["M"].value<double>();
  this->gamma0 = *params["EOS"]["gamma0"].value<double>();
  this->Pc = *params["EOS"]["Pc"].value<double>();
  this->Tc = *params["EOS"]["Tc"].value<double>();

  this->A = 27./64. * this->R0*this->R0 * this->Tc*this->Tc / this->Pc;
  this->B = 1./8. * this->R0 * this->Tc / this->Pc;

  std::string flux = *params["EOS"]["NumFlux"].value<std::string>();
  if (flux == "CENTRAL") {
    this->numerical_flux = CENTRAL;
  } else if (flux == "PEQ_CENTRAL") {
    this->numerical_flux = PEQ_CENTRAL;
  } else {
    std::cout << "Unsupported NumFlux: " << flux << std::endl;
    throw -1;
  }

  this->PrintSelf();

//  this->TestEOS(params);
}

//-------------------------------------------------------------------------------------

void VanDerWaalsEOS::PrintSelf() {
  std::stringstream msg;
  msg << std::endl
      << "------------------------------------" << std::endl
      << "VanDerWaalsEOS::PrintSelf()" << std::endl
      << "------------------------------------" << std::endl;

  msg << "M: " << this->M << std::endl;
  msg << "gamma0: " << this->gamma0 << std::endl;
  msg << "Pc: " << this->Pc << std::endl;
  msg << "Tc: " << this->Tc << std::endl;

  if (this->numerical_flux == CENTRAL)
    msg << "Numerical flux: " << "CENTRAL" << std::endl;
  else if (this->numerical_flux == PEQ_CENTRAL)
    msg << "Numerical flux: " << "PEQ_CENTRAL" << std::endl;

  std::cout << msg.str();
}

//-------------------------------------------------------------------------------------

void VanDerWaalsEOS::CalcNumericalFlux(double R0, double R1,
                                       double RU0, double RU1,
                                       double RE0, double RE1,
                                       double RY0, double RY1,
                                       double P0, double P1,
                                       double T0, double T1,
                                       double &Frho0, double &Frhou0,
                                       double &FrhoE0, double &FrhoY0,
                                       double &Frho1, double &Frhou1,
                                       double &FrhoE1, double &FrhoY1,
                                       double BR0, double BR1,
                                       double BRY0, double BRY1,
                                       double &FA_MARK) {
  double U0 = RU0 / R0;
  double U1 = RU1 / R1;

  double Re0 = RE0 - 0.5*R0*U0*U0;
  double Re1 = RE1 - 0.5*R1*U1*U1;

  if (this->numerical_flux == CENTRAL) {
      double Frho, Frhou, FrhoE, FrhoY;
      double C0f, Mf,PIf, Kf, If, Pf;

      C0f = 0.5*(RU0 + RU1);
      Mf = 0.5*(RU0*U0 + RU1*U1);
      PIf = 0.5*(P0 + P1);
      Kf = 0.5*(RU0*0.5*U0*U0 + RU1*0.5*U1*U1);
      If = 0.5*(Re0*U0 + Re1*U1);
      Pf = 0.5*(U0*P0 + U1*P1);

      Frho = C0f;
      Frhou = Mf + PIf;
      FrhoE = Kf + If + Pf;
      FrhoY = 0.;

      Frho0 = Frho; Frho1 = Frho;
      Frhou0 = Frhou; Frhou1 = Frhou;
      FrhoE0 = FrhoE; FrhoE1 = FrhoE;
      FrhoY0 = FrhoY; FrhoY1 = FrhoY;
  } else if (this->numerical_flux == PEQ_CENTRAL) {
    double Frho, Frhou, FrhoE_, FrhoY;
    double C0f, Mf,PIf, Kf, If0, If1, Pf;

    C0f = 0.5*(RU0 + RU1);
    Mf = 0.5*(RU0*U0 + RU1*U1);
    PIf = 0.5*(P0 + P1);
    Kf = 0.5*(RU0*0.5*U0*U0 + RU1*0.5*U1*U1);
    Pf = 0.5*(U0*P0 + U1*P1);
    // Split flux for energy
    If0 = BR0 * C0f;
    If1 = BR1 * C0f;

    Frho = C0f;
    Frhou = Mf + PIf;
    FrhoE_ = Kf + Pf;
    FrhoY = 0.;

    Frho0 = Frho; Frho1 = Frho;
    Frhou0 = Frhou; Frhou1 = Frhou;
    FrhoE0 = FrhoE_ + If0; FrhoE1 = FrhoE_ + If1;
    FrhoY0 = FrhoY; FrhoY1 = FrhoY;
  } else
    throw -1;
}

//-------------------------------------------------------------------------------------

void VanDerWaalsEOS::RealFluidFromTR() {
  double v = this->M / this->R;

  this->dPdT_v = this->R0 / (v - this->B);
  this->dPdv_T = - this->R0*this->T / (v - this->B)/(v - this->B) + 2.*this->A/v/v/v;

  double cv = this->R0 / (this->gamma0 - 1.);
  double cp = cv - this->T*dPdT_v*dPdT_v/dPdv_T;
  this->gamma = cp / cv;
  this->SoS = std::sqrt(gamma / R / (-1./v/dPdv_T));
  assert(!std::isnan(SoS));

  double dEdR_T = -this->A/this->M/this->M;
  double dEdT_R = this->R0 / (this->gamma0 - 1.) / this->M;
  double dTdR_P = this->M / this->R/this->R * dPdv_T / dPdT_v;
  this->beta_R = this->E + this->R * (dEdR_T + dEdT_R*dTdR_P);

  this->beta_RY = 0.;
}

void VanDerWaalsEOS::SetMixture_TRY(double T, double R, double Y) {
  this->T = T;
  this->R = R;
  this->Y = Y;

  double v = this->M / this->R; // m^3 / kmol
  this->P = this->R0*this->T / (v - this->B) - this->A / v/v;

  double Ehat = this->R0/(this->gamma0 - 1.)*this->T - this->A / v;
  this->E = Ehat / this->M;

  this->RealFluidFromTR();
}

void VanDerWaalsEOS::SetMixture_PRY(double P, double R, double Y) {
  this->P = P;
  this->R = R;
  this->Y = Y;

  double v = this->M / this->R;
  this->T = (this->P + this->A / v/v) * (v - this->B) / this->R0;

  double Ehat = this->R0/(this->gamma0 - 1.)*this->T - this->A / v;
  this->E = Ehat / this->M;

  this->RealFluidFromTR();
}

void VanDerWaalsEOS::SetMixture_ERY(double E, double R, double Y, double Tguess) {
  this->E = E;
  this->R = R;
  this->Y = Y;

  double v = this->M / this->R;
  double Ehat = this->E * this->M;
  this->T = (Ehat + this->A / v) * (this->gamma0 - 1.) / this->R0;

  this->P = this->R0*this->T / (v - this->B) - this->A / v/v;

  this->RealFluidFromTR();
}

void VanDerWaalsEOS::RZFromTP() {
  double a = this->A * this->P / this->R0/this->R0 / this->T/this->T;
  double b = this->B * this->P / this->R0 / this->T;

  double a0 = -a*b;
  double a1 = a;
  double a2 = -(1+b);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; ++i) {
    if (xZ[i] > b)
      Z.push_back(xZ[i]);
  }

  assert(Z.size() == 1);

  double Zo = Z[0];

  double v = Zo * this->R0 * this->T / this->P;
  this->R = this->M / v;
}

void VanDerWaalsEOS::SetMixture_TPY(double T, double P, double Y) {
  this->T = T;
  this->P = P;
  this->Y = Y;

  this->RZFromTP();

  double v = this->M / this->R;
  double Ehat = this->R0/(this->gamma0 - 1.)*this->T - this->A / v;
  this->E = Ehat / this->M;

  this->RealFluidFromTR();
}

void VanDerWaalsEOS::TestEOS(const toml::table& params) {
  double Pin = *params["EOS"]["Pr_in"].value<double>() * this->Pc;
  double Tin = *params["EOS"]["Tr_in"].value<double>() * this->Tc;

  std::cout << "-----------------" << std::endl;
  std::cout << "TP" << std::endl;
  this->SetMixture_TPY(Tin, Pin, 0.);

  std::cout << "T " << this->T << std::endl;
  std::cout << "P " << this->P << std::endl;
  std::cout << "R " << this->R << std::endl;
  std::cout << "E " << this->E << std::endl;
  std::cout << "Sos " << this->SoS << std::endl;
  std::cout << "beta_R " << this->beta_R << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "TR" << std::endl;
  this->SetMixture_TRY(Tin, this->R, 0.);

  std::cout << "T " << this->T << std::endl;
  std::cout << "P " << this->P << std::endl;
  std::cout << "R " << this->R << std::endl;
  std::cout << "E " << this->E << std::endl;
  std::cout << "Sos " << this->SoS << std::endl;
  std::cout << "beta_R " << this->beta_R << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "ER" << std::endl;
  this->SetMixture_ERY(this->E, this->R, 0.);

  std::cout << "T " << this->T << std::endl;
  std::cout << "P " << this->P << std::endl;
  std::cout << "R " << this->R << std::endl;
  std::cout << "E " << this->E << std::endl;
  std::cout << "Sos " << this->SoS << std::endl;
  std::cout << "beta_R " << this->beta_R << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "PR" << std::endl;
  this->SetMixture_PRY(Pin, this->R, 0.);

  std::cout << "T " << this->T << std::endl;
  std::cout << "P " << this->P << std::endl;
  std::cout << "R " << this->R << std::endl;
  std::cout << "E " << this->E << std::endl;
  std::cout << "Sos " << this->SoS << std::endl;
  std::cout << "beta_R " << this->beta_R << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "Step rho" << std::endl;
  double step = 1e-4; double Ro = this->R; double REo = this->R*this->E;
  this->SetMixture_PRY(Pin, this->R*(1.+step), 0.);
  std::cout << "dREdR_P " << (this->R*this->E - REo) / Ro / step << std::endl;;

  throw -1;
}
