#ifndef GENERICEOS_H
#define GENERICEOS_H

#include "toml++/toml.hpp"


class GenericEOS {
public:
  GenericEOS() = default;
  virtual ~GenericEOS();

  virtual void Initialize(const toml::table& params) {
    throw "Not implemented";
  }

  virtual void SetMixture_TPY(double T, double P, double Y) {
    throw "Not implemented";
  }
  virtual void SetMixture_TRY(double T, double R, double Y) {
    throw "Not implemented";
  }
  virtual void SetMixture_ERY(double T, double R, double Y, double Tguess=0.) {
    throw "Not implemented";
  }
  virtual void SetMixture_PRY(double P, double R, double Y) {
    throw "Not implemented";
  }

  virtual void CalcNumericalFlux(double R0, double R1,
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
    throw "Not implemented";
  }

  virtual void PrintSelf();

  virtual double GetP() { return this->P; }
  virtual double GetT() { return this->T; }
  virtual double GetR() { return this->R; }
  virtual double GetE() { return this->E; }
  virtual double GetY() { return this->Y; }
  virtual double GetSoS() { return this->SoS; }
  virtual double GetBeta_R() { return this->beta_R; }
  virtual double GetBeta_RY() { return this->beta_RY; }

  // TODO: other needed variables

protected:
  double E, R, Y;
  double P, T;
  double SoS;

  double beta_R, beta_RY;
};

#endif // GENERICEOS_H
