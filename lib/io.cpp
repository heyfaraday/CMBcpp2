#include <fstream>
#include <complex>
#include <iostream>
#include <vector>

#include "io.hpp"

#include "parameters.hpp"

void o_map(std::string name, long double** map) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            out_file << i << ' ' << j << ' '
                 << std::scientific << map[i][j] << std::endl;
        }
    }

    out_file.close();
}

void i_map(std::string name, long double** map) {

    unsigned int i;
    unsigned int j;
    long double map_;

    std::ifstream in_file;
    in_file.open(name);

    std::string line;

    while(std::getline(in_file, line)) {
        std::stringstream stream(line);
        stream >> i;
        stream >> j;
        stream >> map_;
        map[i][j] = map_;
    }

}

void o_aml(std::string name, long double** aml) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            out_file << m << ' ' << l << ' '
                     << std::scientific << aml[m][l] << std::endl;
        }
    }

    out_file.close();
}

void i_aml(std::string name, long double** aml) {

    unsigned int m;
    unsigned int l;
    long double aml_;

    std::ifstream in_file;
    in_file.open(name);

    std::string line;

    while(std::getline(in_file, line)) {
        std::stringstream stream(line);
        stream >> m;
        stream >> l;
        stream >> aml_;
        aml[m][l] = aml_;
    }
}
