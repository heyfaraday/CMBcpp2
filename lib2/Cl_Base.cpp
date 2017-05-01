#include <cassert>
#include <string>

#include "Cl_Base.hpp"

Cl_Base::Cl_Base(unsigned int n_mod_input) {

    n_mod = n_mod_input;

    array = new long double[n_mod]();
    assert(array != nullptr);
}

Cl_Base::~Cl_Base() {

    delete [] array;
}

unsigned int Cl_Base::n_mod_value() const {
    return n_mod;
}

long double Cl_Base::at(unsigned int l) const {

    assert(l >= 0 and l < n_mod);

    return array[l];
}

void Cl_Base::to(unsigned int l, long double cl) {

    assert(l >= 0 and l < n_mod);

    array[l] = cl;
}

void Cl_Base::fout(std::string name) const {

    FILE* fp = fopen(name.c_str(), "w");

    for (unsigned int l = 0; l < n_mod; ++l) {
        fprintf(fp, "%u %.21Le \n", l, array[l]);
    }

    fclose(fp);
}

void Cl_Base::fin(std::string name) const {

    unsigned int l;
    long double pix;

    FILE* fp = fopen(name.c_str(), "r");

    while(!feof(fp)) {
        fscanf(fp,"%u %Lf \n", &l, &pix);
        assert(l >= 0 and l < n_mod);
        array[l] = pix;
    }

    fclose(fp);
}
