#include <iostream>
#include <fstream> 
#include <chrono>
#include <iomanip> 
#include <cstdio>
#include <mpi.h>
#include <math.h>



#define pi 3.14159265358979323846

constexpr int dx = 1;
constexpr int dz = 1;
constexpr int dt = 1;
constexpr double ro = 1.0;

constexpr int q = 9;           
constexpr int ex[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr int ez[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr double w[]  = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};

/**
 * @brief Parameters for the simulation: for parallelization requirements ONLY CHOOSE nx MULTIPLE OF SIZE.
 * @param dx Grid spacing in x direction
 * @param dz Grid spacing in z direction
 * @param iter Number of iterations
 */
constexpr int nx = 140;      
constexpr int nz = 140;  
constexpr int iter = 50000;


/**
 * @brief Simulation main quantities.
 */

double (*rho)[nz] = (double(*)[nz]) malloc(nx * nz * sizeof(double));
double (*ux)[nz] = (double(*)[nz]) malloc(nx * nz * sizeof(double));
double (*uz)[nz] = (double(*)[nz]) malloc(nx * nz * sizeof(double));
double (*u_abs)[nz] = (double(*)[nz]) malloc(nx * nz * sizeof(double));
double (*u_abs_temp)[nz] = (double(*)[nz]) malloc(nx * nz * sizeof(double));
double (*f)[nx][nz] = (double(*)[nx][nz]) malloc(q * nx * nz * sizeof(double));
double (*f_eq)[nx][nz] = (double(*)[nx][nz]) malloc(q * nx * nz * sizeof(double));
double (*f_temp)[nx][nz] = (double(*)[nx][nz]) malloc(q * nx * nz * sizeof(double));
double (*pressure)[nz] = (double(*)[nz])malloc(nx * nz * sizeof(double));
double (*alpha_old)[nz] = (double(*)[nz])malloc(nx * nz * sizeof(double));