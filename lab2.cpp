#include <mpi.h>

double f(double x) {
    return x;
}
using namespace std;
int main(int argc, char* argv[]) {

    int rank;
    int size;

    double time;
    unsigned long long n = 10000000000;

    double a = 0.0;
    double b = 2.0;

    double sum = 0.0;
    double h = (b - a) / n;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        time = MPI_Wtime();

    for (int i = rank; i < n; i += size)
        sum += f(a + (i + 0.5) * h);
    sum *= h;

    double reduced_sum = 0;
    MPI_Reduce(&sum, &reduced_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        time = MPI_Wtime() - time;
        cout << "Integral value is : " << reduced_sum << endl;
        cout << "Time : " << time << endl;
    }
    MPI_Finalize();
    return 0;
}
