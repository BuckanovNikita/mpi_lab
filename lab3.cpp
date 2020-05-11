#include <mpi.h>
using namespace std;

int get_chunk(int total, int comm_size, int rank) {
    int n = total;
    int q = n / comm_size;
    if (n % comm_size)
        q++;
    int r = comm_size * q - n;
    int chunk = q;
    if (rank >= comm_size - r)
        chunk = q - 1;
    return chunk;
}

int main(int argc, char *argv[]) {
    int n = 4000;
    double t;
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int num_rows = get_chunk(n, comm_size, rank);
    int *rows = static_cast<int *>(malloc(sizeof(int) * num_rows));
    double *a = static_cast<double *>(malloc(sizeof(double) * num_rows * (n + 1)));
    double *x = static_cast<double *>(malloc(sizeof(double) * n));
    double *tmp = static_cast<double *>(malloc(sizeof(double) * (n + 1)));

    for (int i = 0; i < num_rows; i++) {
        rows[i] = rank + comm_size * i;
        srand(rows[i] * (n + 1));
        for (int j = 0; j < n; j++)
            a[i * (n + 1) + j] = rand() % 100 + 1;
        a[i * (n + 1) + n] = rand() % 100 + 1;
    }
    if (rank == 0)
    {
        t = MPI_Wtime();
    }
    int row = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i == rows[row]) {
            MPI_Bcast(&a[row * (n + 1)], n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            for (int j = 0; j <= n; j++)
                tmp[j] = a[row * (n + 1) + j];
            row++;
        } else {
            MPI_Bcast(tmp, n + 1, MPI_DOUBLE, i % comm_size, MPI_COMM_WORLD);
        }
        for (int j = row; j < num_rows; j++) {
            double scaling = a[j * (n + 1) + i] / tmp[i];
            for (int k = i; k < n + 1; k++)
                a[j * (n + 1) + k] -= scaling * tmp[k];
        }
    }
    row = 0;
    for (int i = 0; i < n; i++) {
        x[i] = 0;
        if (i == rows[row]) {
            x[i] = a[row * (n + 1) + n];
            row++;
        }
    }
    row = num_rows - 1;
    for (int i = n - 1; i > 0; i--) {
        if (row >= 0) {
            if (i == rows[row]) {
                x[i] /= a[row * (n + 1) + i];
                MPI_Bcast(&x[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                row--;
            } else
                MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % comm_size, MPI_COMM_WORLD);
        } else
            MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % comm_size, MPI_COMM_WORLD);
        for (int j = 0; j <= row; j++)
            x[rows[j]] -= a[j * (n + 1) + i] * x[i];
    }
    if (rank == 0)
        x[0] /= a[row * (n + 1)];
    MPI_Bcast(x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(tmp);
    free(rows);
    free(a);
    t = MPI_Wtime() - t;
    if (rank == 0) {
        cout << "Time : " << t <<endl;
    }
    free(x);
    MPI_Finalize();
    return 0;
}
