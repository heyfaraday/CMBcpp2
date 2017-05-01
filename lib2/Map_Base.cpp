#include <cassert>
#include <string>

#include "Map_Base.hpp"

Map_Base::Map_Base(unsigned int n_pix_input) {

    n_pix = n_pix_input;

    array = new long double*[n_pix + 1];
    assert(array != nullptr);

    for(unsigned int i = 0; i <= n_pix; ++i) {
        array[i] = new long double[n_pix / 2 + 1]();
        assert(array[i] != nullptr);
    }
}

Map_Base::~Map_Base() {

    for(unsigned int i = 0; i <= n_pix; ++i) {
        delete [] array[i];
    }

    delete [] array;
}

unsigned int Map_Base::n_pix_value() const {
    return n_pix;
}

long double Map_Base::at(unsigned int i, unsigned int j) const {

    assert(i >= 0 and i <= n_pix);
    assert(j >= 0 and j <= n_pix / 2);

    return array[i][j];
}

void Map_Base::to(unsigned int i, unsigned int j, long double pix) {

    assert(i >= 0 and i <= n_pix);
    assert(j >= 0 and j <= n_pix / 2);

    array[i][j] = pix;
}

long double Map_Base::at_border(unsigned int i, unsigned int j) const {

    assert(i >= 0 and i <= n_pix);
    assert(j > 0 and j < n_pix / 2);

    return array[i][j];
}

void Map_Base::to_border(unsigned int i, unsigned int j, long double pix) {

    assert(i >= 0 and i <= n_pix);
    assert(j > 0 and j < n_pix / 2);

    array[i][j] = pix;
}

void Map_Base::fout(std::string name) const {

    FILE* fp = fopen(name.c_str(), "w");

    for (unsigned int i = 0; i <= n_pix; ++i) {
        for (unsigned int j = 0; j <= n_pix / 2; ++j) {
            fprintf(fp, "%u %u %.21Le \n", i, j, array[i][j]);
        }
    }

    fclose(fp);
}

void Map_Base::fin(std::string name) const {

    unsigned int i;
    unsigned int j;
    long double pix;

    FILE* fp = fopen(name.c_str(), "r");

    while(!feof(fp)) {
        fscanf(fp,"%u %u %Lf \n", &i, &j, &pix);
        assert(i >= 0 and i <= n_pix);
        assert(j >= 0 and j <= n_pix / 2);
        array[i][j] = pix;
    }

    fclose(fp);
}
