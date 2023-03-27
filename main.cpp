/*
 Exact Riemann Solver for the Euler Equations

 Assign configuration file with the option -c
   Ex) -c Config/Config.dat

 Developer : Juhyun Kim
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "avocado.h"

// Constants for ideal gas
const double GAMMA = 1.4;
const double GAMMA_P1 = GAMMA + 1;
const double GAMMA_M1 = GAMMA - 1;
const double GAMMA2 = GAMMA * 2;
const double GAMMA_M1P1 = GAMMA_M1 / GAMMA_P1;
const double GAMMA_P1_2 = GAMMA_P1 / GAMMA2;
const double GAMMA_M1_2 = GAMMA_M1 / GAMMA2;
const double GAMMA_inv = 1 / GAMMA;

// Structure to store primitive variables
struct PrimitiveVariables {
  double rho;
  double u;
  double p;
};

// Structure to store conservative variables
struct ConservativeVariables {
  double rho;
  double rho_u;
  double e;
};

// Compute the speed of sound
double sound_speed(double p, double rho) { return std::sqrt(GAMMA * p / rho); }

// Compute the pressure function f(p)
double pressure_function(double p, const PrimitiveVariables& vars) {
  double A = 2 / (GAMMA_P1 * vars.rho);
  double B = GAMMA_M1P1 * vars.p;
  double a = sound_speed(vars.p, vars.rho);

  if (p > vars.p) {
    /// Shock wave
    double sqrt_term = std::sqrt(A / (p + B));
    return (p - vars.p) * sqrt_term;
  } else {
    /// Rarefaction wave
    double power_term = std::pow((p / vars.p), GAMMA_M1_2);
    return 2 * a / GAMMA_M1 * (power_term - 1);
  }
}

// Compute the derivative of the pressure function f'(p)
double pressure_function_prime(double p, const PrimitiveVariables& vars) {
  double A = 2 / (GAMMA_P1 * vars.rho);
  double B = GAMMA_M1P1 * vars.p;
  double a = sound_speed(vars.p, vars.rho);

  if (p > vars.p) {
    /// Shock wave
    double sqrt_term = std::sqrt(A / (p + B));
    return (1 - 0.5 * (p - vars.p) / (B + p)) * sqrt_term;
  } else {
    /// Rarefaction wave
    double power_term = std::pow((p / vars.p), -GAMMA_P1_2);
    return 1 / (vars.p * a) * power_term;
  }
}

// Find the exact solution for the Riemann problem using the Newton-Raphson
// method
double find_exact_pressure(const PrimitiveVariables& left,
                           const PrimitiveVariables& right, double tolerance) {
  double p_guess = 0.5 * (left.p + right.p);
  double max_iterations = 1000;
  double error = 1e6;

  int iter = 0;
  for (iter = 0; iter < max_iterations && error > tolerance; ++iter) {
    double f_left = pressure_function(p_guess, left);
    double f_right = pressure_function(p_guess, right);
    double f = f_left + f_right + right.u - left.u;
    double f_prime = pressure_function_prime(p_guess, left) +
                     pressure_function_prime(p_guess, right);
    double p_new = p_guess - f / f_prime;
    p_new = std::max(tolerance, p_new);

    error = 2.0 * std::abs(p_new - p_guess) / (p_new + p_guess);
    p_guess = p_new;
  }
  if (iter == max_iterations ||
      error > tolerance) {  /// Try again with 1/10 rate if diverged
    //MASTER_MESSAGE("Divergence in Newton-Raphson\n");
    //MASTER_MESSAGE("Trying again with 1/10 step\n");
    for (iter = 0; iter < max_iterations && error > tolerance; ++iter) {
      double f_left = pressure_function(p_guess, left);
      double f_right = pressure_function(p_guess, right);
      double f = f_left + f_right + right.u - left.u;
      double f_prime = pressure_function_prime(p_guess, left) +
                       pressure_function_prime(p_guess, right);
      double p_new = p_guess - 0.1 * f / f_prime;
      p_new = std::max(tolerance, p_new);

      error = 2.0 * std::abs(p_new - p_guess) / (p_new + p_guess);
      p_guess = p_new;
    }
    if (iter == max_iterations || error > tolerance) {
      ERROR_MESSAGE("Divergence in Newton-Raphson\n");
    }
  }
  return p_guess;
}

// Compute the exact solution for the Riemann problem
void solve_riemann(const PrimitiveVariables& left,
                   const PrimitiveVariables& right, double t,
                   PrimitiveVariables& result, double x, double tolerance) {
  double p_star = find_exact_pressure(left, right, tolerance);
  double u_star =
      0.5 * (left.u + right.u) + 0.5 * (pressure_function(p_star, right) -
                                        pressure_function(p_star, left));
  if (x / t < u_star) {
    /// Left side of the contact discontinuity
    if (p_star > left.p) {
      /// Shock wave
      double shock_speed =
          left.u - sound_speed(left.p, left.rho) *
                       std::sqrt(GAMMA_P1_2 * p_star / left.p + GAMMA_M1_2);
      if (x / t < shock_speed) {
        result = left;
      } else {
        result.rho = left.rho * (p_star / left.p + GAMMA_M1P1) /
                     (GAMMA_M1P1 * (p_star / left.p) + 1);
        result.u = u_star;
        result.p = p_star;
      }
    } else {
      /// Rarefaction wave
      const double a = sound_speed(left.p, left.rho);
      double head_speed = left.u - a;
      double tail_speed = u_star - a * std::pow(p_star / left.p, GAMMA_M1_2);
      if (x / t < head_speed) {
        result = left;
      } else if (x / t > tail_speed) {
        result.rho = left.rho * std::pow((p_star / left.p), GAMMA_inv);
        result.u = u_star;
        result.p = p_star;
      } else {
        const double a = sound_speed(left.p, left.rho);
        const double temp = 2 / GAMMA_P1 + GAMMA_M1P1 / a * (left.u - x / t);
        result.rho = left.rho * std::pow(temp, 2 / GAMMA_M1);
        result.u = 2 / GAMMA_P1 * (a + 0.5 * GAMMA_M1 * left.u + x / t);
        result.p = left.p * std::pow(temp, 2 * GAMMA / GAMMA_M1);
      }
    }
  } else {
    /// Right side of the contact discontinuity
    if (p_star > right.p) {
      /// Shock wave
      double shock_speed =
          right.u + sound_speed(right.p, right.rho) *
                        std::sqrt((GAMMA_P1_2 * p_star / right.p + GAMMA_M1_2));
      if (x / t > shock_speed) {
        result = right;
      } else {
        result.rho = right.rho * (p_star / right.p + GAMMA_M1P1) /
                     (GAMMA_M1P1 * (p_star / right.p) + 1);
        result.u = u_star;
        result.p = p_star;
      }
    } else {
      /// Rarefaction wave
      const double a = sound_speed(right.p, right.rho);
      double head_speed = right.u + a;
      double tail_speed = u_star + a * std::pow(p_star / right.p, GAMMA_M1_2);

      if (x / t > head_speed) {
        result = right;
      } else if (x / t < tail_speed) {
        result.rho = right.rho * std::pow((p_star / right.p), GAMMA_inv);
        result.u = u_star;
        result.p = p_star;
      } else {
        const double a = sound_speed(right.p, right.rho);
        const double temp = 2 / GAMMA_P1 - GAMMA_M1P1 / a * (right.u - x / t);
        result.rho = right.rho * std::pow(temp, 2 / GAMMA_M1);
        result.u = 2 / GAMMA_P1 * (-a + 0.5 * GAMMA_M1 * right.u + x / t);
        result.p = right.p * std::pow(temp, 2 * GAMMA / GAMMA_M1);
      }
    }
  }
}

int main(int argc, char* argv[]) {
  std::vector<std::string> add_argv;
  avocado::AVOCADO_Initialize(argc, argv, "Exact Riemann Solver", add_argv);

  // Configuration
  const auto& config = AVOCADO_CONFIG;
  const double domain_length = std::stod(
      config->GetConfigValue("DomainLength"));  /// domain length from x=0
  const double discont_pos =
      std::stod(config->GetConfigValue("DiscontinuityPosition"));
  const int num_cells = std::stoi(config->GetConfigValue("NumCells"));
  const double gamma = std::stod(config->GetConfigValue("Gamma"));
  const double target_time =
      std::stod(config->GetConfigValue("Time"));  /// Target time
  const double d_L = std::stod(config->GetConfigValue("DensityL"));
  const double u_L = std::stod(config->GetConfigValue("VelocityL"));
  const double p_L = std::stod(config->GetConfigValue("PressureL"));
  const double d_R = std::stod(config->GetConfigValue("DensityR"));
  const double u_R = std::stod(config->GetConfigValue("VelocityR"));
  const double p_R = std::stod(config->GetConfigValue("PressureR"));
  const double tolerance = std::stod(config->GetConfigValue("Tolerance"));

  PrimitiveVariables left = {
      d_L, u_L, p_L};  /// Left initial conditions (density, velocity, pressure)
  PrimitiveVariables right = {
      d_R, u_R,
      p_R};  /// Right initial conditions (density, velocity, pressure)

  std::ofstream outfile;
  std::string filename = config->GetConfigValue("ReturnDir");
  avocado::MakeDirectory(filename);
  filename += config->GetConfigValue("ExportFile");
  outfile.open(filename, std::ios::trunc);
  if (!outfile.is_open()) ERROR_MESSAGE(filename + " is not opened!");
  outfile << "TITLE = \"Exact Riemann Solution \"" << std::endl;
  outfile << "VARIABLES = \"x\",\"d\",\"u\",\"p\"" << std::endl;

  const double dx = domain_length / num_cells;
  std::vector<double> x_coords(num_cells);
  PrimitiveVariables result;
  const double start_x = 0.5 * dx - discont_pos;
  for (int icell = 0; icell < num_cells; icell++) {
    const double x = start_x + dx * icell;

    solve_riemann(left, right, target_time, result, x, tolerance);
    outfile << x + discont_pos << "\t" << result.rho << "\t" << result.u << "\t"
            << result.p << "\t" << std::endl;
  }
  outfile << std::endl;
  outfile.close();

  MASTER_MESSAGE("Program normal exit.\n");
  avocado::AVOCADO_Finalize();

  return 0;
}