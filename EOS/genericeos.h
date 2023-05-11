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
  virtual void SetMixture_ERY(double T, double R, double Y) {
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
                                 double& Frho, double& Frhou,
                                 double& FrhoE, double& FrhoY) {
    throw "Not implemented";
  }

  virtual void PrintSelf();

  virtual double GetP() { return this->P; }
  virtual double GetT() { return this->T; }
  virtual double GetR() { return this->R; }
  virtual double GetE() { return this->E; }
  virtual double GetY() { return this->Y; }
  virtual double GetSoS() { return this->SoS; }

  // TODO: other needed variables

protected:
  double E, R, Y;
  double P, T;
  double SoS;
};

#endif // GENERICEOS_H
