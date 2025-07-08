

using namespace problem_setting;

void write_to_file(const std::string& filename, const std::string& content) {
    std::ofstream file(filename, std::ios::app); // Open file in append mode
    if (file.is_open()) {
        file << content << std::endl;
        file.close();
    } else {
        std::cerr << "Failed to open " << filename << " for writing." << std::endl;
    }
}

void exchange_halos(double (*f_temp)[nx][nz], int q, int local_nx, int rank, int size) {
    MPI_Status status;

    for (int i = 0; i < q; i++) {
        if(rank < size - 1)
            MPI_Send(&f_temp[i][(rank + 1)*local_nx - 1][0], nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        if(rank > 0)
            MPI_Recv(&f_temp[i][(rank - 1)*local_nx + local_nx - 1][0], nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        if(rank > 0)
            MPI_Send(&f_temp[i][rank*local_nx][0], nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        if(rank < size - 1)
            MPI_Recv(&f_temp[i][(rank + 1)*local_nx][0], nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
        if(rank == 0){
            MPI_Send(&f_temp[i][0][0], nz, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD);
        }
        if(rank == size - 1){
            MPI_Recv(&f_temp[i][0][0], nz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        if(rank == size - 1){
            MPI_Send(&f_temp[i][nx - 1][0], nz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        if(rank == 0){
            MPI_Recv(&f_temp[i][nx - 1][0], nz, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, &status);
        }
    }
}


