#include "pengrobinsoneos.h"
#include "utils.h"


#include <iostream>
#include <algorithm>
#include <pugixml.hpp>

//--------------------------------------------------------------------------------

void PengRobinsonEOS::Initialize(const toml::table &params) {
  this->gas_constant_universal = 8314.4621; // [J kmol^-1 K^-1]
  this->N_Av = 6.022e26; // kmol^-1
  this->Boltzman = 1.380649e-23;  // J/K

  const auto& stringList = *params["EOS"]["Species"].as_array();
  for (const auto& element : stringList) {
    this->species.push_back(*element.value<std::string>());
  }
  assert(this->species.size() == n_species);

  // Allocate arrays
  this->Xvec.resize(this->n_species);
  this->Yvec.resize(this->n_species);
  this->MW.resize(this->n_species);

  this->Tcrit.resize(this->n_species);
  this->Pcrit.resize(this->n_species);
  this->rhocrit.resize(this->n_species);
  this->Vcrit.resize(this->n_species);
  this->Zcrit.resize(this->n_species);
  this->omega.resize(this->n_species);

  this->NasaCoef = 7;
  this->nasa_poly_coeff.resize(this->NasaCoef * this->n_species * 2);
  this->nasa_poly_bounds.resize(this->n_species * 2);
  this->hspecies_ig.resize(this->n_species);
  this->cpspecies_ig.resize(this->n_species);
  this->sspecies_ig.resize(this->n_species);

  this->Tcrit_IJ.resize(this->n_species * this->n_species);
  this->Vcrit_IJ.resize(this->n_species * this->n_species);
  this->Zcrit_IJ.resize(this->n_species * this->n_species);
  this->Pcrit_IJ.resize(this->n_species * this->n_species);
  this->omega_IJ.resize(this->n_species * this->n_species);

  this->cst_b.resize(this->n_species);
  this->cst_a.resize(this->n_species * this->n_species);
  this->cst_c.resize(this->n_species * this->n_species);

  this->A_IJ.resize(this->n_species * this->n_species);
  this->dAmdN.resize(this->n_species);
  this->d2AmdTdN.resize(this->n_species);
  this->dK1dN.resize(this->n_species);
  this->dPdN.resize(this->n_species);
  this->dVdN.resize(this->n_species);

  this->hspecies.resize(this->n_species);
  this->sspecies.resize(this->n_species);
  this->muspecies.resize(this->n_species);
  this->edspecies.resize(this->n_species);


  this->ReadNasaPolynomials();
  this->ReadCriticalProperties();
  this->SetRealFluidConstants();
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::ReadNasaPolynomials() {
  std::string nasa_poly_dbpath = NASAPOLYPATH;

  std::vector<double> readFromXMLFile_low;
  std::vector<double> readFromXMLFile_high;

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(nasa_poly_dbpath.c_str());

  LOOP_l_N(this->n_species) {
    // search for species
    pugi::xml_node
        species = xmlDoc.child("ctml").child("speciesData").first_child();
    while (species.attribute("name").value() != this->species[l]) {
      species = species.next_sibling();
    }

    pugi::xml_node lowTemp = species.child("thermo").first_child();
    pugi::xml_node highTemp = lowTemp.next_sibling();

    // read temperature ranges
    this->nasa_poly_bounds[0 + 2 * l].first =
        lowTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[0 + 2 * l].second =
        highTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[1 + 2 * l].first =
        lowTemp.attribute("Tmax").as_double();
    this->nasa_poly_bounds[1 + 2 * l].second =
        highTemp.attribute("Tmax").as_double();

    readFromXMLFile_low =
        Tokenize(lowTemp.child("floatArray").child_value(), " ,\n");
    readFromXMLFile_high =
        Tokenize(highTemp.child("floatArray").child_value(), " ,\n");

    LOOP_k_N(this->NasaCoef) {
      this->nasa_poly_coeff[k + this->NasaCoef * l] = readFromXMLFile_low[k];
      this->nasa_poly_coeff[k + this->NasaCoef * l
          + this->NasaCoef * this->n_species] = readFromXMLFile_high[k];
    }
  }
}

void PengRobinsonEOS::ReadCriticalProperties() {
  LOOP_l_N(this->n_species) {
    if (this->species[l] == "O2") {
      this->MW[l] = 31.9988; //O2
      this->Tcrit[l] = 154.5800;
      this->Pcrit[l] = 5.0430e+6;
      this->rhocrit[l] = 436.140;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.0222;
    } else if (this->species[l] == "H2") {
      this->MW[l] = 2.01588; //H2
      this->Tcrit[l] = 33.1450;
      this->Pcrit[l] = 1.2964e+6;
      this->rhocrit[l] = 31.262;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = -0.219;
    } else if (this->species[l] == "H2O") {
      this->MW[l] = 18.01528;
      this->Tcrit[l] = 647.096;
      this->Pcrit[l] = 22.0640e+6;
      this->rhocrit[l] = 322.0;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.3443;
    } else if (this->species[l] == "OH") {
      this->MW[l] = 17.0073; //OH kg/kmol
      this->Tcrit[l] = 0.0;
      this->Pcrit[l] = 0.0;
      this->rhocrit[l] = 0.0;
      this->Vcrit[l] = 0.0;
      this->Zcrit[l] = 0.0;
      this->omega[l] = 0.0;
    } else if (this->species[l] == "N2") {
      this->MW[l] = 28.0134; //N2 kg/kmol
      this->Tcrit[l] = 126.1900;
      this->Pcrit[l] = 3.3958e+6;
      this->rhocrit[l] = 313.3;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.03720;
    } else if (this->species[l] == "NH3") {
      this->MW[l] = 17.031; //NH3 kg/kmol
      this->Tcrit[l] = 405.40;
      this->Pcrit[l] = 11.3330e+6;
      this->rhocrit[l] = 225.000;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.25601;
    } else if (this->species[l] == "C8H18,isooctane") {
      this->MW[l] = 114.23;
      this->Tcrit[l] = 543.9;
      this->Pcrit[l] = 25.7e5;
      this->rhocrit[l] = 244.4522;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.394;
    } else if (this->species[l] == "NC12H26") {
      this->MW[l] = 170.33484;
      this->Tcrit[l] = 658.10;
      this->Pcrit[l] = 1.817e+6;
      this->rhocrit[l] = 226.5453372;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.574;
    } else if (this->species[l] == "CH4") {
      this->MW[l] = 16.04;
      this->Tcrit[l] = 190.6;
      this->Pcrit[l] = 46.1e5;
      this->rhocrit[l] = 162.0;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.011;
    } else if (this->species[l] == "C7H16,n-heptane") {
      this->MW[l] = 100.2019;
      this->Tcrit[l] = 540;
      this->Pcrit[l] = 27.4e5;
      this->rhocrit[l] = 2.35*100.2019;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.349;
    } else {
      std::cout << " WARNING -> Unknown species :[" << this->species[l]
           << "]. No critical properties found." << std::endl;
      throw -1;
    }
  }
}

void PengRobinsonEOS::SetRealFluidConstants() {
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      double tmp_k = 0; // zero binary-interaction parameter
      // if (k != l)
      //   tmp_k = 0.156;

      int apos = k * this->n_species + l;

      this->Tcrit_IJ[apos] = std::sqrt(this->Tcrit[l] * this->Tcrit[k]) * (1.0 - tmp_k);
      this->Vcrit_IJ[apos] = std::pow(std::pow(this->Vcrit[l], 1.0 / 3.0) + std::pow(this->Vcrit[k], 1.0 / 3.0), 3.0) / 8.0;
      this->Zcrit_IJ[apos] = 0.5 * (this->Zcrit[l] + this->Zcrit[k]);
      this->Pcrit_IJ[apos] = this->Zcrit_IJ[apos] * this->gas_constant_universal * this->Tcrit_IJ[apos] / this->Vcrit_IJ[apos];
      this->omega_IJ[apos] = 0.5 * (this->omega[l] + this->omega[k]);
    }
  }

  LOOP_k_N(this->n_species) {
    this->cst_b[k] = 0.077796 * this->gas_constant_universal * this->Tcrit[k] / this->Pcrit[k];

    LOOP_l_N(this->n_species) {
      int apos = k * n_species + l;
      this->cst_a[apos] = 0.457236 * std::pow(this->gas_constant_universal * this->Tcrit_IJ[apos], 2.0)  / this->Pcrit_IJ[apos];
      if (this->omega_IJ[apos] <=0.49)
        this->cst_c[apos] = 0.37464 + 1.54226 * this->omega_IJ[apos] - 0.26992 * std::pow(this->omega_IJ[apos], 2);
      else
        this->cst_c[apos] = 0.379642 + 1.485030*this->omega_IJ[apos] - 0.164423*std::pow(this->omega_IJ[apos], 2) + 0.016666*std::pow(this->omega_IJ[apos], 3);
    }
  }
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncRealFluidThermodynamicsFromTemperatureDensity() {
  double v = this->MW_M / this->R;

  this->Am = 0.;
  this->Bm = 0.;
  this->dAmdT = 0.;
  this->d2AmdT2 = 0.;

  LOOP_k_N(this->n_species) {
    this->dAmdN[k] = 0.;
    this->d2AmdTdN[k] = 0.;

    this->Bm += this->Xvec[k] * this->cst_b[k];

    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->Xvec[l] * this->Xvec[k];
      this->A_IJ[apos] = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])), 2);
      double G = this->cst_c[apos] * std::sqrt(this->T / this->Tcrit_IJ[apos])
          / (1.0 + cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])));
      double D = this->cst_c[apos] * (1.0 + cst_c[apos]) * this->Tcrit_IJ[apos]
          / this->Pcrit_IJ[apos] * std::sqrt(this->Tcrit_IJ[apos] / this->T);
      this->Am += X_X * this->A_IJ[apos];
      this->dAmdT -= X_X * this->A_IJ[apos] * G;
      this->d2AmdT2 += X_X * D;

      this->dAmdN[k] += this->Xvec[l] * this->A_IJ[apos];
      this->d2AmdTdN[k] += this->Xvec[l] * this->A_IJ[apos] * G;
    }
    this->dAmdN[k] *= 2.0;
    this->d2AmdTdN[k] *= -2.0 / this->T;
  }
  this->dAmdT /= this->T;
  this->d2AmdT2 *= 0.457236 * std::pow(this->gas_constant_universal, 2) / (2.0 * this->T);
  this->dPdT = this->gas_constant_universal / (v - this->Bm)
      - this->dAmdT / (std::pow(v, 2) + 2.0 * v * this->Bm - std::pow(this->Bm, 2));
  double arg = this->gas_constant_universal * this->T * (v + this->Bm)
      * std::pow((v / (v - this->Bm) + this->Bm / (v + this->Bm)), 2);
  this->dPdV = -this->gas_constant_universal * this->T / std::pow((v - this->Bm), 2)
      * (1.0 - 2.0 * this->Am / arg);
  this->expansivity = -this->dPdT / (v * this->dPdV); //ideal gas: equal to 3.34E-3 (1/K)
  this->K1 = 1.0 / (std::sqrt(8.0) * this->Bm)
      * std::log((v + (1 - std::sqrt(2.0)) * this->Bm) / (v + (1 + std::sqrt(2.0)) * this->Bm));

  double temp = v * v + 2.0 * this->Bm * v - this->Bm * this->Bm;
  LOOP_k_N(this->n_species) {
    this->dPdN[k] = this->gas_constant_universal * this->T / (v - this->Bm) +
        this->gas_constant_universal * this->T * this->cst_b[k]
            / std::pow((v - this->Bm), 2) - this->dAmdN[k] / temp
        + 2.0 * this->Am * this->cst_b[k] * (v - this->Bm) / std::pow(temp, 2);
    this->dVdN[k] = -this->dPdN[k] / this->dPdV;
    this->dK1dN[k] = 1.0 / temp * this->dVdN[k] - this->cst_b[k] / this->Bm *
        (this->K1 + v / temp);
  }
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncIdealFluidThermodynamicsFromTemperature() {
  double T1 = this->T;
  double T2 = this->T * this->T;
  double T3 = T2 * this->T;
  double T4 = T3 * this->T;

  this->h_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4 / 5.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
          + this->NasaCoef * this->n_species] / T1;
    }
    this->h_ig += this->Xvec[l] * (hspecies_ig[l]);
  }
  this->h_ig *= (this->T * this->gas_constant_universal) / this->MW_M;

  this->cp_ig = 0.;
  this->cv_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T <= 1000.0) {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4;
    } else {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species];
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4;
    }
    double tmp = cpspecies_ig[l] * this->gas_constant_universal;
    this->cp_ig += this->Xvec[l] * tmp; //cp_mol
    this->cv_ig += this->Xvec[l] * (tmp - this->gas_constant_universal);
  }
  this->cp_ig /= this->MW_M;
  this->cv_ig /= this->MW_M;

  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0] * std::log(T1);
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 2.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 3.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 4.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 6];
    } else {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species] * std::log(T);
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2 / 2.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3 / 3.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4 / 4.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 6
          + this->NasaCoef * this->n_species];
    }
    this->sspecies_ig[l] -= std::log(this->Xvec[l]*this->P/1e5);
    if (this->Xvec[l] < 1e-10)
      sspecies_ig[l] = 0.;
    if (std::isnan(sspecies_ig[l])) {
      std::cout << "Something wrong Ideal. T = " << this->T << ", rho = " << R
                << ", YF = " << this->Yvec[0] << " , P = " << this->P << std::endl;
      throw -1;
    }
  }
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncRZFromPressureTemperature() {
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  if (Z.size() == 1) {
    this->Z = Z[0];
    this->R = this->P / this->gas_constant / this->T / this->Z;
    return;
  } else {
    std::vector<double> lnPhi(Z.size());

    int iMin = 0;
    for (int i = 0; i < Z.size(); i++) {
      lnPhi[i] = -std::log(Z[i] - B)
          - A / B / std::sqrt(8) * std::log((Z[i] + (1 + std::sqrt(2)) * B) / (Z[i] + (1 - std::sqrt(2)) * B))
          + Z[i] - 1;
      if (lnPhi[i] < lnPhi[iMin]) iMin = i;
    }
    this->Z =  Z[iMin];
    this->R = this->P / this->gas_constant / this->T / this->Z;
    return;
  }
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncPZFromTemperatureDensity() {
  double v = this->MW_M / this->R; // m^3 / kmol

  this->P = (this->gas_constant_universal * this->T) / (v - this->Bm)
      - this->Am / (std::pow(v, 2) + 2.0 * v * this->Bm - pow(this->Bm, 2));

  this->Z = this->P * v / this->gas_constant_universal / this->T;

  double Tcmix, Pcmix, wcmix, Zcmix;
  this->ComputeMixtureCrit(Tcmix, Pcmix, wcmix, Zcmix);

  if (this->P < 0 || this->CheckThreeRoots()) {
    // Psat estimate
    double Psat = Pcmix * std::pow(10.0, (7. / 3. * (1. + wcmix)) * (1. - Tcmix / this->T));

    double a = this->Am;
    double b = this->Bm;

    double A = a * Psat / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
    double B = b * Psat / this->gas_constant_universal / this->T;

    double a0 = -(A * B - B * B - B * B * B);
    double a1 = A - 3 * B * B - 2 * B;
    double a2 = -(1 - B);

    std::vector<double> Z(3);
    int n = SolveP3(&Z[0], a2, a1, a0);
    if (n == 1)
      this->P = Psat; // short-cut estimate wrong with PR
    else {
      std::sort(Z.begin(), Z.end());

      double rhoL = Psat / this->gas_constant / this->T / Z[0];
      double rhoV = Psat / this->gas_constant / this->T / Z[2];

      if (this->P < 0 || (this->R < rhoL && this->R > rhoV))
        this->P = Psat;
    }
  }
  if (this->P < 0 || std::isnan(this->P)) {
    std::cout << "PZ: Impossible pressure for diffusion: " << this->P << ", " << this->P << std::endl;
    throw -1;
  }
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncEFromTemperatureDensity() {
  double dep = (this->Am - this->T * this->dAmdT) * this->K1 / this->MW_M;
  this->E = this->h_ig - this->gas_constant * this->T + dep;

  double departureCv = (-this->T * this->d2AmdT2 * this->K1) / this->MW_M;
  double departureCp = departureCv
      + (-this->T * std::pow(this->dPdT, 2) / this->dPdV - this->gas_constant_universal)
          / this->MW_M;
  this->cv = this->cv_ig + departureCv;
  this->cp = this->cp_ig + departureCp;
  this->gamma = this->cp / this->cv;
  double v = this->MW_M / this->R;
  double isocompressibility = -1. / (v * this->dPdV);
  this->SoS = std::sqrt(this->gamma / (this->R * std::fabs(isocompressibility)));

}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SyncPartialPropertiesFromTemperatureDensity() {
  // TODO: Don't think this is needed?
  // hspecies_ig is molar, without rt
  double temp = this->Am - this->T * this->dAmdT;
  double rt = this->gas_constant_universal * this->T;

  LOOP_k_N(this->n_species) {
    this->edspecies[k] = this->dK1dN[k] * temp +
        this->K1 * (this->dAmdN[k] - this->T * this->d2AmdTdN[k]);
    this->hspecies[k] = this->hspecies_ig[k] * rt - rt + this->edspecies[k] +
        this->P * dVdN[k];
    this->hspecies[k] /= this->MW[k];
    this->edspecies[k] /= this->MW[k];
  }

  // muspecies, sspecies
  LOOP_k_N(this->n_species) {
    double logphik = this->cst_b[k] / this->Bm * (this->Z - 1.) - std::log(this->Z - this->Bm*this->P/rt)
                              - this->Am / this->Bm / rt / 2.8284 * std::log((this->Z+2.4142*this->Bm*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm);
    this->muspecies[k] = this->hspecies_ig[k]*rt - this->sspecies_ig[k]*rt + rt * logphik;
    this->muspecies[k] /= this->MW[k];
    this->sspecies[k] = (this->hspecies[k] - this->muspecies[k]) / this->T;
  }
//  if (std::isnan(muspecies[0])) {
//    std::cout << "Something wrong Partial. T = " << this->T << ", rho = " << this->rho
//              << ", YF = " << this->Y[0] << " , P = " << this->P << ", " << this->ZForDiff - this->Bm*this->PForDiff/rt << ", " << (this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt) << ", " << temp1
//              << ", " << this->cst_b[0] / this->Bm * (this->Z - 1.) - std::log(this->ZForDiff - this->Bm*this->PForDiff/rt)
//              << ", " << this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[0]/this->Am - this->cst_b[0]/this->Bm)
//              << ", " << this->cst_b[0] / this->Bm * (this->Z - 1.) - std::log(this->ZForDiff - this->Bm*this->PForDiff/rt)
//                 - this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[0]/this->Am - this->cst_b[0]/this->Bm) << std::endl;
//    throw -1;
//  }
}

//-------------------------------------------------------------------------------------

bool PengRobinsonEOS::CheckThreeRoots() {
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  return Z.size() > 1;
}

void PengRobinsonEOS::ComputeMixtureCrit(double &Tcmix, double &Pcmix, double &wcmix, double &Zcmix) {
  Tcmix = 0.;
  Pcmix = 0.;
  wcmix = 0.;
  Zcmix = 0.;
  double Vcmix = 0.;

  LOOP_k_N(this->n_species) {
    Tcmix += this->Xvec[k] * this->Tcrit[k];
    Zcmix += this->Xvec[k] * this->Zcrit[k];
    Vcmix += this->Xvec[k] * this->Vcrit[k];
    wcmix += this->Xvec[k] * this->omega[k];
  }
  Pcmix = Zcmix * this->gas_constant_universal * Tcmix / Vcmix;
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::SetMixture_TRY(double T, double R, double Y) {
  this->T = T;
  this->R = R;
  this->Y = Y;

  this->SetMassFractionFromY();

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
}

void PengRobinsonEOS::SetMixture_TPY(double T, double P, double Y) {
  this->T = T;
  this->P = P;
  this->Y = Y;

  this->SetMassFractionFromY();

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncRZFromPressureTemperature();

  this->SetMixture_TRY(T, R, Y);
}

void PengRobinsonEOS::SetMixture_ERY(double E, double R, double Y, double Tguess) {
  this->E = E;
  this->R = R;
  this->Y = Y;

  this->SetMassFractionFromY();

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->PrepareRealFluidThermodynamicsForTFromE();
  this->T = GetTemperatureFromEnergy(E, Tguess);

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
}

//-------------------------------------------------------------------------------------

void PengRobinsonEOS::PrepareRealFluidThermodynamicsForTFromE() {
  double v = this->MW_M / this->R;

  this->Bm = 0.;

  LOOP_k_N(this->n_species) {
    this->Bm += this->Xvec[k] * this->cst_b[k];
  }
  this->K1 = 1.0 / (std::sqrt(8.0) * this->Bm)
      * std::log((v + (1 - std::sqrt(2.0)) * this->Bm) / (v + (1 + std::sqrt(2.0)) * this->Bm));
}

double PengRobinsonEOS::GetTemperatureFromEnergy(const double E_in, double T_guess) {
  double T_0 = 300.;
  if (T_guess > 50. && T_guess < 4000.)
    T_0 = T_guess;
  double E_0 = this->EFromTFunc(T_0);

  double T_1 = (E_0 > E_in) ? 0.999 * T_0 : 1.001 * T_0;
  double E_1 = this->EFromTFunc(T_1);

  int counter = 0;
//  while (std::fabs(E_1 - E_in) > 1e-8) {
  while (std::fabs(T_1 - T_0) > 1e-10) {
    double T_tmp = T_1;
    double E_tmp = E_1;

    T_1 -= (E_1 - E_in) * (T_1 - T_0) / (E_1 - E_0);
    T_1 = std::max(T_1, 20.0);
    E_1 = this->EFromTFunc(T_1);

    T_0 = T_tmp;
    E_0 = E_tmp;
    counter++;
    if (counter > 1000) {
      std::cout << "Over 1000 iter" << std::endl;
      return T_guess;
    }
  }

  return T_1;
}

double PengRobinsonEOS::EFromTFunc(double T_in) {
  // Update Pengrobinson Factors
  this->Am = 0.;
  this->dAmdT = 0.;
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->Xvec[l] * this->Xvec[k];
      this->A_IJ[apos] = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(T_in / this->Tcrit_IJ[apos])), 2);
      double G = this->cst_c[apos] * std::sqrt(T_in / this->Tcrit_IJ[apos])
          / (1.0 + cst_c[apos] * (1.0 - std::sqrt(T_in / this->Tcrit_IJ[apos])));
      this->Am += X_X * this->A_IJ[apos];
      this->dAmdT -= X_X * this->A_IJ[apos] * G;
    }
  }
  this->dAmdT /= T_in;

  // Find Ideal and real energy contribution
  double T1 = T_in;
  double T2 = T_in * T_in;
  double T3 = T2 * T_in;
  double T4 = T3 * T_in;

  this->h_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (T_in < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4 / 5.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
          + this->NasaCoef * this->n_species] / T1;
    }
    this->h_ig += this->Xvec[l] * (hspecies_ig[l]);
  }
  this->h_ig *= (T_in * this->gas_constant_universal) / this->MW_M;

  double dep = (this->Am - T_in * this->dAmdT) * this->K1 / this->MW_M;
  this->E = this->h_ig - this->gas_constant * T_in + dep;

  return this->E;
}






