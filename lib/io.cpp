#include <fstream>

#include "io.hpp"

#include "parameters.hpp"

void io_map(std::string name, long double** map) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix/2; ++j) {
            out_file << i << ' ' << j << ' '
                 << std::scientific << map[i][j] << std::endl;
        }
    }

    out_file.close();
}
