#include "macros.hpp"

namespace problem_setting
{

constexpr const char *base_path = "results/";

// Adimensionalize the problem

constexpr double L = 0.15; // Characteristic length
constexpr double C_t = 0.01; // Characteristic time
constexpr double C_l = L / nx; // Characteristic length in adimensional units
constexpr double v_ref = 6.66e-3;
constexpr double V = v_ref * C_t / C_l; // Reference velocity in adimensional units
constexpr double nu = 1e-6;
constexpr double Re = v_ref * L / nu; // Reynolds number
//constexpr double tau = 0.5 + 3*(nx * V) / Re;
constexpr double nu_adim = nu * C_t / (C_l * C_l); 
constexpr double beta = 1 / (6 * nu_adim + 1);


double const ux_ref(double x, double z)
{
    return 0;
}

double const uz_ref(double x, double z)
{
    return 0;
}

}






