#include <iostream>
#include <cmath>
#include <fstream>

#include <fft.hpp>
#include <utils.hpp>

#include <parameters.hpp>
#include <io.hpp>

#include <fftw3.h>

int main() {

    std::cout << "Hello, World!" << std::endl;

    typedef std::numeric_limits<double> dbl;

    double d = 3.14159265358979412;
    double a = 1.0;
    int b = 100;
    std::cout.precision(dbl::max_digits10);
    std::cout << "Pi: " << std::scientific << d << std::endl;
    std::cout << "a: " << std::scientific << sin(a * d / 2.0) << std::endl;
    std::cout << "int b: " << std::scientific << b << std::endl;

    std::ofstream myfile;
    myfile.open("out.dat");
    myfile << "Writing this to a file.\n";
    myfile.close();

    fftwl_complex *in, *out;
    fftwl_plan p;
    fftwl_plan p1;

    int N = 8;

    in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(fftwl_complex) * static_cast<unsigned long int>(N)));
    out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(fftwl_complex) * static_cast<unsigned long int>(N)));

//    in2 = (fftw*) fftwl_malloc((fftwl) * N);

    p = fftwl_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//
//    fftw_plan p3 = fftw_plan_dft_r2c_1d()

    for (int i = 0; i < N; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    in[0][0] = 0.0L;
    in[1][1] = 1.0L;
//    in[1][1] = 0.0;
//    in[2][0] = 0.0;
//    in[3][0] = 0.0;
//    in[4][0] = 0.0;

    for (int i = 0; i < N; ++i) {
        std::cout << "i=" << i << ": " << in[i][1] << std::endl;
    }

    std::cout << std::endl;

    fftwl_execute(p);

    for (int i = 0; i < N; ++i) {
        std::cout << "i=" << i << ": " << out[i][0] << std::endl;
        out[i][1] = 0.0L;
    }
    std :: cout << std::endl;


    fftwl_destroy_plan(p);

    p1 = fftwl_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftwl_execute(p1);

    for (int i = 0; i < N; ++i) {
        std::cout << "i=" << i << ": " << in[i][1] / 8.0L << std::endl;
    }

    fftwl_destroy_plan(p1);

    int num = 0;

    for (int j = 0; j < N; ++j) {
        if (out[j][0] > 0.0L)
        {
            num = num + 1;
        }
    }

    std::cout << num << std::endl;

    long double** matrix = n_matrix_generator(10, 10);

    for (int i = 0; i < 10; ++i) {
        std::cout << matrix[1][i] << std::endl;
    }

    n_matrix_destroyer(matrix, 10);

    std::cout << "nmod: " << nmod << std::endl;

    return 0;
}
