#include "Solvers/genericsolver.h"
#include "toml++/toml.hpp"

int main() {
  toml::table params = toml::parse_file("input.toml");
  GenericSolver solver;
  solver.Initialize(params);
  solver.RunSolver();
  return 0;
}
