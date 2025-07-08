#include "parameters.hpp"
#include "helper_functions.hpp"
#include "kernel.hpp"

int main(int argc, char *argv[]) {

    using namespace problem_setting;

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int local_nx = nx / size;    

    for(int i = 0; i < nx; i++){
        for(int k = 0; k < nz; k++){
            ux[i][k] = ux_ref(i + 0.5, k + 0.5);
            uz[i][k] = uz_ref(i + 0.5, k + 0.5);
            rho[k][i] = 1.0;
            alpha_old[k][i] = 2.0; 
        }        
    }

    for (int i = 0 + rank * local_nx; i < local_nx * (rank + 1); i++) {
        for (int k = 0; k < nz; k++) {
            for (int a = 0; a < q; a++) {
                double f_x = (2 - pow(1.0 + 3.0 * ux[i][k] * ux[i][k], 0.5)) * pow(((2.0 * ux[i][k] + pow(1.0 + 3.0 * ux[i][k] * ux[i][k], 0.5)) / (1 - ux[i][k])), ex[a]);
                double f_z = (2 - pow(1.0 + 3.0 * uz[i][k] * uz[i][k], 0.5)) * pow(((2.0 * uz[i][k] + pow(1.0 + 3.0 * uz[i][k] * uz[i][k], 0.5)) / (1 - uz[i][k])), ez[a]);
                f_eq[a][i][k]   = w[a]*rho[i][k] * f_x * f_z;
                f_temp[a][i][k] = f_eq[a][i][k];
            }
        }        
    }

    for (int i = 0; i < q; i++) {
        MPI_Allgather(&f_temp[i][rank * local_nx][0], local_nx * nz, MPI_DOUBLE, &f_temp[i][0][0], local_nx * nz, MPI_DOUBLE, MPI_COMM_WORLD);
    } 

    // Optional: ANSI escape codes for colored text

    if(rank == 0) {
    
    const std::string cyan = "\033[36m";
    const std::string reset = "\033[0m";
    std::cout << cyan << "=== Parameters ===" << reset << "\n";
    const std::string red = "\033[31m";
    // Print variable names in cyan, values in red
    std::cout << std::setw(25) << "Characteristic length (L): " << red << L << reset << " m\n";
    std::cout << std::setw(25) << "Characteristic time interval(C_t): " << red << C_t << reset << " s\n";
    std::cout << std::setw(25) << "Characteristic space step(C_t): " << red << C_l << reset << " m\n";
    std::cout << std::setw(25) << "Grid size (nx): " << red << nx << reset << "\n";
    std::cout << std::setw(25) << "Grid size (nz): " << red << nz << reset << "\n";
    std::cout << std::setw(25) << "Total grid points: " << red << nx * nz << reset << "\n";
    std::cout << std::setw(25) << "Reference velocity (v_ref): " << red << v_ref << reset << " m/s\n";
    std::cout << std::setw(25) << "Kinematic viscosity (nu): " << red << nu << reset << "\n";
    std::cout << std::setw(25) << "Reynolds number (Re): " << red << Re << reset << "\n";
    std::cout << cyan << "======================================" << reset << "\n";
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Starting simulation" << std::endl;

  
    for(int n=0; n<iter; n++){
    stream(rank, local_nx);
    compute_macroscopic_quantities(rank, local_nx);
    compute_equilibrium(rank, local_nx);
    calculate_alpha(rank, local_nx);
    collide(rank, local_nx);
    
    // Exchange halos
    if(size > 1){
        exchange_halos(f_temp, q, local_nx, rank, size);
    }

    compute_convergence(rank, local_nx, n);

    if(n % 100 == 0 || n == iter - 1) {

        char filename[128];


        // Allgather for ux
        MPI_Allgather(&ux[rank * local_nx][0], local_nx * nz, MPI_DOUBLE, &ux[0][0], local_nx * nz, MPI_DOUBLE, MPI_COMM_WORLD);

        // Allgather for uz
        MPI_Allgather(&uz[rank * local_nx][0], local_nx * nz, MPI_DOUBLE, &uz[0][0], local_nx * nz, MPI_DOUBLE, MPI_COMM_WORLD);

        // Allgather for u_abs
        MPI_Allgather(&u_abs[rank * local_nx][0], local_nx * nz, MPI_DOUBLE, &u_abs[0][0], local_nx * nz, MPI_DOUBLE, MPI_COMM_WORLD);

        sprintf(filename, "results/ux_%d.vtr", n);
        write_vtr_scalar(filename, "ux", ux, nx, nz, dx, dz);

        sprintf(filename, "results/uz_%d.vtr", n);
        write_vtr_scalar(filename, "uz", uz, nx, nz, dx, dz);

        sprintf(filename, "results/u_abs_%d.vtr", n);
        write_vtr_scalar(filename, "u_abs", u_abs, nx, nz, dx, dz);
    }


    }



    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Simulation completed in " << elapsed_time << " seconds." << std::endl;
    MPI_Finalize();
   
    return 0;
}

