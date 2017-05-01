#include <iostream>
//
#include "Map_Base.hpp"
#include "Alm_Base.hpp"
#include "Cl_Base.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

//    Map_Base* map = new Map_Base(8);
    Map_Base map(8);
    Alm_Base alm(30);

    std::cout << map.n_pix_value() << std::endl;
    std::cout << map.at(8, 4) << std::endl;

    std::cout << alm.n_mod_value() << std::endl;
    std::cout << alm.at_cos(0, 0) << std::endl;

    map.fin("test_fout.dat");

    map.fout("test2_fout.dat");

    alm.fin("test3_fout.dat", "test4_fout.dat");

    Cl_Base cl(8);

    cl.fin("test4.dat");
    std::cout << cl.at(7) << std::endl;

    map.to(3, 4, 1203);

//    delete map;

    return 0;
}
