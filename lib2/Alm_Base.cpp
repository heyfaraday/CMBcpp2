#include <cassert>
#include <string>

#include "Alm_Base.hpp"

Alm_Base::Alm_Base(unsigned int n_mod_input) {

    n_mod = n_mod_input;

    array_cos = new long double*[n_mod];
    assert(array_cos != nullptr);

    array_sin = new long double*[n_mod];
    assert(array_sin != nullptr);

    for(unsigned int l = 0; l < n_mod; ++l) {
        array_cos[l] = new long double[l + 1]();
        assert(array_cos[l] != nullptr);
    }

    for(unsigned int l = 0; l < n_mod; ++l) {
        array_sin[l] = new long double[l + 1]();
        assert(array_sin[l] != nullptr);
    }
}

Alm_Base::~Alm_Base() {

    for(unsigned int l = 0; l < n_mod; ++l) {
        delete [] array_cos[l];
    }

    delete [] array_cos;

    for(unsigned int l = 0; l < n_mod; ++l) {
        delete [] array_sin[l];
    }

    delete [] array_sin;
}

unsigned int Alm_Base::n_mod_value() const {
    return n_mod;
}

long double Alm_Base::at_cos(unsigned int l, unsigned int m) const {

    assert(l >= 0 and l < n_mod);
    assert(m >= 0 and m <= l);

    return array_cos[l][m];
}

void Alm_Base::to_cos(unsigned int l, unsigned int m, long double alm) {

    assert(l >= 0 and l < n_mod);
    assert(m >= 0 and m <= l);

    array_cos[l][m] = alm;
}

long double Alm_Base::at_sin(unsigned int l, unsigned int m) const {

    assert(l >= 0 and l < n_mod);
    assert(m >= 0 and m <= l);

    return array_sin[l][m];
}

void Alm_Base::to_sin(unsigned int l, unsigned int m, long double alm) {

    assert(l >= 0 and l < n_mod);
    assert(m >= 0 and m <= l);

    array_sin[l][m] = alm;
}

void Alm_Base::fout(std::string name_cos, std::string name_sin) const {

    FILE* fp_cos = fopen(name_cos.c_str(), "w");
    FILE* fp_sin = fopen(name_sin.c_str(), "w");

    for (unsigned int l = 0; l < n_mod; ++l) {
        for (unsigned int m = 0; m <= l; ++m) {
            fprintf(fp_cos, "%u %u %.21Le \n", l, m, array_cos[l][m]);
        }
    }

    for (unsigned int l = 0; l < n_mod; ++l) {
        for (unsigned int m = 0; m <= l; ++m) {
            fprintf(fp_sin, "%u %u %.21Le \n", l, m, array_sin[l][m]);
        }
    }

    fclose(fp_cos);
    fclose(fp_sin);
}

void Alm_Base::fin(std::string name_cos, std::string name_sin) const {

    unsigned int l;
    unsigned int m;
    long double pix;

    FILE* fp_cos = fopen(name_cos.c_str(), "r");
    FILE* fp_sin = fopen(name_sin.c_str(), "r");

    while(!feof(fp_cos)) {
        fscanf(fp_cos,"%u %u %Lf \n", &l, &m, &pix);
        assert(l >= 0 and l < n_mod);
        assert(m >= 0 and m <= l);
        array_cos[l][m] = pix;
    }

    while(!feof(fp_sin)) {
        fscanf(fp_sin,"%u %u %Lf \n", &l, &m, &pix);
        assert(l >= 0 and l < n_mod);
        assert(m >= 0 and m <= l);
        array_sin[l][m] = pix;
    }

    fclose(fp_cos);
    fclose(fp_sin);
}
