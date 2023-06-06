#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <iostream>

#include "toml++/toml.hpp"

// EOS
#include "EOS/genericeos.h"
#include "EOS/idealgaseos.h"
#include "EOS/vanderwaalseos.h"

static GenericEOS* MakeEOS(const toml::table& params) {
  GenericEOS* ret;

  std::string eos_type = *params["EOS"]["type"].value<std::string>();
  if (eos_type == "IdealGas") {
    IdealGasEOS* pt = new IdealGasEOS;
    pt->Initialize(params);
    ret = pt;
  } else if (eos_type == "VanDerWaals") {
    VanDerWaalsEOS* pt = new VanDerWaalsEOS;
    pt->Initialize(params);
    ret = pt;
  } else {
    std::cout << "Unsupported EOS: " << eos_type << std::endl;
    throw -1;
  }

  return ret;
}

#endif // INTERFACE_H
