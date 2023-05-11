#ifndef IDEALGASEOS_H
#define IDEALGASEOS_H

#include "genericeos.h"

class IdealGasEOS : public GenericEOS {
public:
  enum NumFlux {CENTRAL, PEQ_CENTRAL};
  IdealGasEOS() = default;
  ~IdealGasEOS() override = default;

  virtual void Initialize(const toml::table& params) override;

  virtual void SetMixture_TPY(double T, double P, double Y) override;
  virtual void SetMixture_TRY(double T, double R, double Y) override;
  // E is just specific internal energy
  virtual void SetMixture_ERY(double E, double R, double Y) override;
  virtual void SetMixture_PRY(double P, double R, double Y) override;

  virtual void CalcNumericalFlux(double R0, double R1,
                                 double RU0, double RU1,
                                 double RE0, double RE1,
                                 double RY0, double RY1,
                                 double P0, double P1,
                                 double T0, double T1,
                                 double& Frho, double& Frhou,
                                 double& FrhoE, double& FrhoY) override;

  virtual void PrintSelf() override;

private:
  virtual void CalcMixingRule();

  NumFlux numerical_flux;

  double gamma0, gamma1;
  double M0, M1;
  constexpr static const double R0 = 8314; // J/kmol/K

  double M_bar, R_bar;
  double gamma_bar;
};

#endif // IDEALGASEOS_H
