#include <iostream>
#include <omp.h>
#include <ctime>
#include <cmath>

using namespace std;

const double L = 1; // sizes of the plate

double f(double x, double y) {
    return 2.0 * L * (y + x) - 2 * (x * x + y * y);
}

double u_(double x, double y) {
    return x * y * (L - x) * (L - y);
}

double CalculateError(int n, double *approx) {
    double h = L / n;
    double max_err = 0;
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++) {
            double cur_err = fabs(approx[i * (n - 1) + j] - u_((i + 1) * h, (j + 1) * h));
            if (cur_err > max_err)
                max_err = cur_err;
        }
    return max_err;
}

void SerialFFT1D(complex<double> *inputSignal, complex<double> *outputSignal,
                 int size, int step, int offset) {


}

void bitReversing()

void four1(double *data, unsigned long nn) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    // reverse-binary reindexing
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            swap(data[j - 1], data[i - 1]);
            swap(data[j], data[i]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };

    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = -(2 * M_PI / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j - 1] - wi * data[j];
                tempi = wr * data[j] + wi * data[j - 1];

                data[j - 1] = data[i - 1] - tempr;
                data[j] = data[i] - tempi;
                data[i - 1] += tempr;
                data[i] += tempi;
            }
            wtemp = wr;
            wr += wr * wpr - wi * wpi;
            wi += wi * wpr + wtemp * wpi;
        }
        mmax = istep;
    }
}

void ComputeDecision(double *approx, int n, void (*SinFT)(double *, double *, int), int numThreads) {

    double h = L / n;
    omp_set_num_threads(numThreads);
    double **f_four = new double *[n - 1];
    for (int i = 0; i < n - 1; i++) {
        f_four[i] = new double[n - 1];
    }

#pragma omp parallel for
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1; j++) {
            approx[i * (n - 1) + j] = f((i + 1) * h, (j + 1) * h)
            //* mult1 * h;
        }
        (*SinFT)(&approx[i * (n - 1)], f_four[i], n);
        for (int j = 0; j < n - 1; j++) {
            f_four[i][j] = -f_four[i][j];
        }
    }
// transpose matrix

#pragma omp parallel for
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1)

    }
    double **dl = new double *[numThreads];
    double **d = new double *[numThreads];
    double **du = new double *[numThreads];

    for (int i = 0; i < numThreads; i++) {
        dl[i] = new double[n - 2];
        d[i] = new double[n - 1];
        du[i] = new double[n - 2];
    }

    double sudiag_val = 1.0 / (h * h);
    for (int k = 0; k < n - 2; k++) {
        for (int i = 0; i < n - 2; i++) {
            dl[i] = sudiag_val;
            du[i] = sudiag_val;
        }
        double tmp = sin(M_PI * (k + 1) / (2 * n));
        double l_k = 4 / (h * h) * tmp * tmp;
        double diag_val = -l_k - 2.0 / (h * h);

        for (int i = 0; i < n - 1; i++)
            d[i] = diag_val;
    }
}


int main() {
    omp_set_num_threads(4);
    int num_threads = 4;
    int n = 10;

    double *approx = new double[(n - 1) * (n - 1)];
    clock_t start = clock();
    ComputeDecision(approx, n, FFT, num_threads);
    clock_t finish = clock();
    time = (double) (finish - start) / CLOCKS_PER_SEC;
    double err = CalculateError(n, approx);
    cout << n << time << err << endl;
    delete[] approx;

    return 0;
}