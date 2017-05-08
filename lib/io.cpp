#include <fstream>
#include <complex>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "io.hpp"

#include "parameters.hpp"

void o_map(std::string name, long double** map) {

//    std::ofstream out_file;
//    out_file.open(name);
//    typedef std::numeric_limits<long double> dbl;
//    out_file.precision(dbl::max_digits10);

    FILE* fp = fopen(name.c_str(), "w");

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            fprintf(fp, "%u %u %.21Le \n", i, j, map[i][j]);
        }
    }
    fclose(fp);
}

void o_map_norm(std::string name, long double** map) {

    FILE* fp = fopen(name.c_str(), "w");

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            fprintf(fp, "%.21Le %.21Le %.21Le \n", i * map_parameter, j * map_parameter, map[i][j]);
        }
    }
    fclose(fp);
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
            out_file << m << "  " << l << "  "
                     << std::scientific << aml[m][l] << std::endl;
        }
    }

    out_file.close();
}

void i_aml(std::string name, long double** aml) {

    unsigned int m;
    unsigned int l;
    long double aml_;

//    std::ifstream in_file;
//    in_file.open(name);
//    std::string line;

    FILE* fp = fopen(name.c_str(), "r");

    while(!feof(fp)) {
        fscanf(fp,"%u %u %Lf \n", &m, &l, &aml_);
        if (m < nmod and l < nmod) {
            aml[m][l] = aml_;
        }
    }

    fclose(fp);
}

void o_cl(std::string name, long double* cl) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    for (unsigned int l = 0; l < nmod; ++l) {
        out_file << l << "  "
                 << std::scientific << cl[l] << std::endl;
    }

    out_file.close();

}

void i_cl(std::string name, long double* cl) {

    unsigned int l;
    long double cl_;

    std::ifstream in_file;
    in_file.open(name);

    std::string line;

    while(std::getline(in_file, line)) {
        std::stringstream stream(line);
        stream >> l;
        stream >> cl_;
        if (l < nmod) {
            cl[l] = cl_;
        }
    }
}
