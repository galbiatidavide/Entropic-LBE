// Stream with periodic boundary conditions using a linkwise approach (non optimized maintain for safety)

void stream(int rank, int local_nx){

    int b  = 0;   
    int ix, iz;
            

    for (int i = 0 + rank*local_nx; i < local_nx*(rank+1); i++){
        for (int k=0; k<nz; k++) {                
            for (int a=0; a<q; a++) {

                if(a == 0)
                    b = 0; 
                else if (a == 1)
                    b = 3;
                else if (a == 2)
                    b = 4;
                else if (a == 3)
                    b = 1;
                else if (a == 4)
                    b = 2;
                else if (a == 5)
                    b = 7;
                else if (a == 6)
                    b = 8;
                else if (a == 7)
                    b = 5;
                else if (a == 8)
                    b = 6;

                ix = i-ex[a]; 
                iz = k-ez[a];


                if(ix < 0 || ix > nx-1 || iz < 0 || iz > nz-1) {
                    if(iz < 0){
                        double vx;
                        if(i < 70)
                            vx = V - V*exp(-0.1 * i);
                        if(i >= 70)
                            vx = V - V*exp(-0.1 * (nx - 1 - i));

                        f[a][i][k] = f_temp[b][i][k] - 6*w[b]*rho[i][k]*(ex[b]*vx);
                    }
                    else 
                        f[a][i][k] = f_temp[b][i][k];
                }
                else 
                    f[a][i][k] = f_temp[a][ix][iz];
            }                    
        }        
    }
}

void compute_macroscopic_quantities(int rank, int local_nx){
    for (int i = 0 + rank*local_nx; i < local_nx*(rank+1); i++) {
        for (int k=0; k<nz; k++) {
            rho[i][k] = 0;
            ux[i][k]  = 0;
            uz[i][k]  = 0;
            
            for (int a=0; a<q; a++) {
                rho[i][k] =  rho[i][k] + f[a][i][k];
                ux[i][k]  =  ux[i][k]  + ex[a]*f[a][i][k];
                uz[i][k]  =  uz[i][k]  + ez[a]*f[a][i][k];
            }
            u_abs_temp[i][k] = u_abs[i][k];
            ux[i][k] = (ux[i][k])/rho[i][k];
            uz[i][k] = (uz[i][k])/rho[i][k];
            u_abs[i][k] = sqrt(ux[i][k]*ux[i][k] + uz[i][k]*uz[i][k]);
        }
    }
}

void compute_equilibrium(int rank, int local_nx){
    for (int i=0 + rank*local_nx; i<local_nx*(rank+1); i++) {
        for(int k=0; k<nz; k++) {
            for (int a =0 ; a<q; a++) {
                double f_x = (2 - pow(1.0 + 3.0 * ux[i][k] * ux[i][k], 0.5)) * pow(((2.0 * ux[i][k] + pow(1.0 + 3.0 * ux[i][k] * ux[i][k], 0.5)) / (1 - ux[i][k])), ex[a]);
                double f_z = (2 - pow(1.0 + 3.0 * uz[i][k] * uz[i][k], 0.5)) * pow(((2.0 * uz[i][k] + pow(1.0 + 3.0 * uz[i][k] * uz[i][k], 0.5)) / (1 - uz[i][k])), ez[a]);
                f_eq[a][i][k]   = w[a]*rho[i][k] * f_x * f_z;
            }
        }        
    }
}


void collide(int rank, int local_nx) {

    for(int i = 0 + rank*local_nx; i < local_nx*(rank+1); i++){ 
        for (int k = 0; k < nz; k++) {
            for (int a = 0; a < q; a++) {
                double omega = -(f[a][i][k]- f_eq[a][i][k])*alpha_old[i][k]*beta;
                f_temp[a][i][k] = f[a][i][k]+omega;
            }
        }        
    }
    
}

void compute_convergence(int rank, int local_nx, int n) {
    double temp1_local = 0.0;
    double temp2_local = 0.0;

    for(int i = 0 + rank*local_nx; i < local_nx*(rank+1); i++){ 
        for (int k = 0; k < nz; k++) {
            double u_curr = u_abs[i][k];
            double u_prev = u_abs_temp[i][k];
            double diff = u_curr - u_prev;

            temp1_local += diff * diff;
            temp2_local += u_curr * u_curr;
        }
        
    }

    double temp1_global = 0.0, temp2_global = 0.0;
    MPI_Reduce(&temp1_local, &temp1_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&temp2_local, &temp2_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double err_con = std::sqrt(temp1_global / temp2_global);
        std::cout << "[# Iter " << n << "] Convergence error is " << err_con << std::endl;
    }
}

void entropy_function(int x, int y, double &H, double &H_der, double &H_alpha, double &H_der_alpha) {
        
    for (int a = 0; a<q; a++) {
        H     += f[a][x][y]*log(f[a][x][y]/w[a]);
        H_der += 1+log(f[a][x][y]/w[a]);
    }         
    for (int a = 0; a<q; a++) {
        double f_new = f[a][x][y]*(1-alpha_old[x][y]) + alpha_old[x][y]*f_eq[a][x][y];
        H_alpha     += f_new*log(f_new/w[a]);
        H_der_alpha += (1-alpha_old[x][y])*(1+log(f_new/w[a]));
    }         
}

void calculate_alpha(int rank, int local_nx){
         
    double alpha_new = 0.0;
    double error = 1.0;
         
    do {
            
        for(int i = 0 + rank*local_nx; i < local_nx*(rank+1); i++){ 
            for (int k = 0; k < nz; k++) {
                    
                double H = 0;
                double H_der = 0;
                double H_alpha = 0;
                double H_der_alpha = 0;
                entropy_function(i,k, H, H_der, H_alpha, H_der_alpha);
                double H_function     = H - H_alpha;
                double H_function_der = H_der - H_der_alpha;
                alpha_new = alpha_old[i][k] - H_function/H_function_der;
                if (alpha_new == 0.0) {
                    double temp_array[q];
                    for (int a=0; a<q; a++) {
                        temp_array[a] = abs(f[a][i][k]/(f_eq[a][i][k]-f[a][i][k]));
                    }
                    double max_alpha = temp_array[0];
                    for (int a=1; a<q; a++) {
                        if (temp_array[a]>temp_array[a-1]) {
                            max_alpha = temp_array[a];
                        }
                    }
                    alpha_new = max_alpha;
                }                
                
                error = abs(alpha_new-alpha_old[i][k]); 
                alpha_old[i][k] = alpha_new; 
            }
        }
            
    } while (error>1e-6);
         
}

void write_vtr_scalar(const char *filename, const char *varname,
                      double (*field)[nz], int nx, int nz, double dx, double dz) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    fprintf(fp,
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        "<RectilinearGrid WholeExtent=\"0 %d 0 0 0 %d\">\n"
        "<Piece Extent=\"0 %d 0 0 0 %d\">\n"
        "<PointData Scalars=\"%s\">\n"
        "<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",
        nx-1, nz-1, nx-1, nz-1, varname, varname
    );

    // Data values
    for (int j = 0; j < nz; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fp, "%f ", field[i][j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp,
        "</DataArray>\n"
        "</PointData>\n"
        "<Coordinates>\n"
        "<DataArray type=\"Float64\" Name=\"XCoordinates\" format=\"ascii\">\n"
    );
    for (int i = 0; i < nx; i++) fprintf(fp, "%f ", i * dx);
    fprintf(fp, "\n</DataArray>\n");

    fprintf(fp,
        "<DataArray type=\"Float64\" Name=\"YCoordinates\" format=\"ascii\">\n"
        "0.0\n</DataArray>\n");

    fprintf(fp,
        "<DataArray type=\"Float64\" Name=\"ZCoordinates\" format=\"ascii\">\n");
    for (int j = 0; j < nz; j++) fprintf(fp, "%f ", j * dz);
    fprintf(fp,
        "\n</DataArray>\n"
        "</Coordinates>\n"
        "</Piece>\n"
        "</RectilinearGrid>\n"
        "</VTKFile>\n"
    );

    fclose(fp);
}







































