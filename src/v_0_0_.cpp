#include <iostream>
#include <math.h>
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

    fftw_complex *in, *out;
    fftw_plan p;

    int N = 8;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
        out[i][0] = 0.0;
        out[i][1] = 0.0;
    }

    // in[0][0] = 1.0 / sqrt(4.0 * d);
    in[1][1] = 1.0;

    fftw_execute(p);

    for (int i = 0; i < N; ++i) {
        std::cout << "i=" << i << ": " << out[i][0] << std::endl;
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    int num = 0;

    for (int j = 0; j < N; ++j) {
        if (out[j][0] > 0.0)
        {
            num = num + 1;
        }
    }

    std::cout << num << std::endl;

    double** matrix = n_matrix_generator(10, 10);

    for (int i = 0; i < 10; ++i) {
        std::cout << matrix[1][i] << std::endl;
    }

    n_matrix_destroyer(matrix, 10);

    std::cout << "nmod: " << nmod << std::endl;

    return 0;
}
