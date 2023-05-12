#ifndef PENGROBINSONEOS_H
#define PENGROBINSONEOS_H

#include "genericeos.h"

#define LOOP_l_N(n) for(int l = 0; l < n; ++l)
#define LOOP_k_N(n) for(int k = 0; k < n; ++k)

class PengRobinsonEOS : public GenericEOS {
public:
  PengRobinsonEOS() = default;
  ~PengRobinsonEOS() = default;

  virtual void Initialize(const toml::table& params) override;

  virtual void SetMixture_TPY(double T, double P, double Y) override;
  virtual void SetMixture_TRY(double T, double R, double Y) override;
  virtual void SetMixture_ERY(double E, double R, double Y, double Tguess) override;
//  virtual void SetMixture_PRY(double P, double R, double Y) override;

  virtual void CalcNumericalFlux(double R0, double R1,
                                 double RU0, double RU1,
                                 double RE0, double RE1,
                                 double RY0, double RY1,
                                 double P0, double P1,
                                 double T0, double T1,
                                 double& Frho, double& Frhou,
                                 double& FrhoE, double& FrhoY,
                                 double BR0, double BR1,
                                 double BRY0, double BRY1,
                                 double& FA_MARK) override;

  virtual void PrintSelf() override;

private:
  inline void SetMolecularWeightMixtureFromY() {
    this->MW_M = 0.;
    LOOP_l_N(this->n_species)
        this->MW_M += this->Yvec[l] / MW[l];
    this->MW_M = 1. / this->MW_M;

    this->gas_constant = this->gas_constant_universal / this->MW_M;
  }

  inline void SetMolarFractionFromY() {
    LOOP_l_N(this->n_species)
        this->Xvec[l] = this->Yvec[l] * this->MW_M / this->MW[l];
  }

  inline void SetMassFractionFromY() {
    this->Yvec[0] = this->Y;
    this->Yvec[1] = 1. - this->Y;
  }


  void ReadNasaPolynomials();
  void ReadCriticalProperties();
  void SetRealFluidConstants();

  void SyncRealFluidThermodynamicsFromTemperatureDensity();
  void SyncIdealFluidThermodynamicsFromTemperature();

  void SyncRZFromPressureTemperature();
  void SyncPZFromTemperatureDensity();
  void SyncEFromTemperatureDensity();
  void SyncPartialPropertiesFromTemperatureDensity();

  // For iterating to find T from E
  void PrepareRealFluidThermodynamicsForTFromE();
  double GetTemperatureFromEnergy(const double E_in, double Tguess);
  double EFromTFunc(double T);

  // Utilities
  // Test if current T and P can give more than one solution
  bool CheckThreeRoots();
  void ComputeMixtureCrit(double& Tcmix, double& Pcmix, double& wcmix, double& Zcmix);


  double gas_constant_universal, N_Av, Boltzman;

  constexpr static const int n_species = 2;
  std::vector<std::string> species;

  std::vector<double> Xvec;
  std::vector<double> Yvec;

  double MW_M;
  double gas_constant;

  double Z;

  double cp, cv, gamma;

  std::vector<double> A_IJ;
  double Am;
  double Bm;
  double dAmdT;
  double d2AmdT2;
  std::vector<double> dAmdN;
  std::vector<double> d2AmdTdN;
  double K1;
  std::vector<double> dK1dN;

//  double BT;
//  double betaT;
  double expansivity;
  double dPdT;
  double dPdV;
  std::vector<double> dPdN;
  std::vector<double> dVdN;

  std::vector<double> hspecies;
  std::vector<double> sspecies;
  std::vector<double> muspecies;
  std::vector<double> edspecies;


  std::vector<double> MW;

  std::vector<double> Tcrit;
  std::vector<double> Pcrit;
  std::vector<double> rhocrit;
  std::vector<double> Vcrit;
  std::vector<double> Zcrit;
  std::vector<double> omega;

  int NasaCoef;
  std::vector<double> nasa_poly_coeff;
  std::vector<std::pair<double, double>> nasa_poly_bounds;
  std::vector<double> cpspecies_ig;
  std::vector<double> hspecies_ig;
  std::vector<double> sspecies_ig;
  double h_ig;
  double cp_ig;
  double cv_ig;

  std::vector<double> Tcrit_IJ;
  std::vector<double> Vcrit_IJ;
  std::vector<double> Zcrit_IJ;
  std::vector<double> Pcrit_IJ;
  std::vector<double> omega_IJ;

  std::vector<double> cst_a;
  std::vector<double> cst_b;
  std::vector<double> cst_c; // k
};

#endif // PENGROBINSONEOS_H
