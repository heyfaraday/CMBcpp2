#include <fstream>

#include "io.hpp"

#include "parameters.hpp"

void io_map(std::ofstream file, double** map) {

    typedef std::numeric_limits<double> dbl;

    file.precision(dbl::max_digits10);

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix/2; ++j) {
            file << i << ' ' << j << ' '
                 << std::scientific << map[i][j];
        }
    }
}
