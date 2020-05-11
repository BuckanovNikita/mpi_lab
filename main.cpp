#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t n;
    if (argc >= 2)
        n = strtoul(argv[1], nullptr, 10);
    else
    {
        cout <<"Need n as input arg"<<endl;
        return 1;
    }

    double *a, *c;
    double *b = (double*)malloc(n * n * sizeof(double));
    double *working_a = (double*)malloc(n * n * sizeof(double));
    double *working_c = (double*)malloc(n * n * sizeof(double));



    size_t num_workers = size;
    if (num_workers > n)
        num_workers = n;

    size_t block_size = num_workers > 0 ? n / num_workers : 0;
    size_t mod_block_size = n - num_workers * block_size;

    double start_time;
    if (rank == 0) {
        a = (double*)malloc(n * n * sizeof(double));
        c = (double*)malloc(n * n * sizeof(double));
        ifstream in_a("a.txt");
        ifstream in_b("b.txt");
        for (int i = 0; i < n * n; i++) {
            in_a >> a[i];
            in_b >> b[i];
        }
        start_time = MPI_Wtime();
    }
    MPI_Bcast(b, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a, block_size * n, MPI_DOUBLE, working_a, block_size * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int row = 0; row < block_size; row++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += working_a[row * n + j] * b[j * n + i];
            working_c[row * n + i] = sum;
        }
    }

    MPI_Gather(working_c, block_size * n, MPI_DOUBLE, c, block_si
    ze*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        if (mod_block_size > 0) {
            int start_row = block_size * num_workers;
            for (int row = start_row; row < n; row++) {
                for (int i = 0; i < n; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < n; j++)
                        sum += a[row * n + j] * b[j * n + i];
                    c[row * n + i] = sum;
                }
            }
        }

        double total_time = MPI_Wtime() - start_time;
        cout << "Matrix size: " << n <<endl;
        cout << "Total time: " << total_time <<" sec." <<endl;
    }

    MPI_Finalize();
    if (rank == 0)
    {
        ifstream in_c("c.txt");
        double err = 0;
        for(int i=0; i<n*n; i++)
        {
            double true_c;
            in_c >> true_c;
            err  = max(err, abs(true_c - c[i]));
        }
        cout << "Max error: " << err <<endl;
    }
    return 0;
}
