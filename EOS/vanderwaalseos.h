#ifndef VANDERWAALSEOS_H
#define VANDERWAALSEOS_H

#include "genericeos.h"


class VanDerWaalsEOS : public GenericEOS {
 public:
  enum NumFlux { CENTRAL, PEQ_CENTRAL };
  VanDerWaalsEOS() = default;
  ~VanDerWaalsEOS() override = default;

  virtual void Initialize(const toml::table& params) override;

  virtual void SetMixture_TPY(double T, double P, double Y) override;
  virtual void SetMixture_TRY(double T, double R, double Y) override;
  // E is just specific internal energy
  virtual void SetMixture_ERY(double E, double R, double Y, double Tguess=0.) override;
  virtual void SetMixture_PRY(double P, double R, double Y) override;

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
                                 double& FA_MARK) override;

  virtual void PrintSelf() override;

 private:
  void TestEOS(const toml::table& params);

  void RealFluidFromTR();
  void RZFromTP();

  constexpr static const double R0 = 8314; // J/kmol/K

  double M, gamma0;
  double Pc, Tc;
  double A, B;

  double dPdT_v, dPdv_T;
  double gamma;

  NumFlux numerical_flux;
};

#endif // VANDERWAALSEOS_H
