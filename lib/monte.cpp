#include <string>
#include <cmath>
#include <fstream>
#include <iostream>

#include "functionals.hpp"
#include "parameters.hpp"
#include "constants.hpp"

#include "monte.hpp"

void circles_area(long double** map, long double level, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double theta;
    unsigned int j_center = 0;
    long double v1, v2;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        v1 = area(map, level, npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1);
        v2 = area(map, level, npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        out_file << std::scientific << v1 << " "
                 << npix / 4 + j_center << " "
                 << npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1
                 << std::endl
                 << std::scientific << v2 << " "
                 << npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)) << " "
                 << npix / 4 - j_center - 1
                 << std::endl;

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }

    out_file.close();
}

void circles_length(long double** map, long double level, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double theta;
    unsigned int j_center = 0;
    long double v1, v2;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        v1 = length(map, level, npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1);
        v2 = length(map, level, npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        out_file << std::scientific << v1 << " "
                 << npix / 4 + j_center << " "
                 << npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1
                 << std::endl
                 << std::scientific << v2 << " "
                 << npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)) << " "
                 << npix / 4 - j_center - 1
                 << std::endl;

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }

    out_file.close();
}

void circles_points(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                    std::string name) {

    long double theta;
    unsigned int j_center = 0;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        points_classifier(map_x, map_y, cos_lm, sin_lm, name, npix / 4 + j_center, npix / 4 + j_center +
                static_cast<unsigned int>(hring * sinl(theta)) - 1);

        points_classifier(map_x, map_y, cos_lm, sin_lm, name, npix / 4 -
                j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }
}

void circles_singular_points(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                             long double** u_cos_lm, long double** u_sin_lm, std::string name) {

    long double theta;
    unsigned int j_center = 0;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        singular_points_classifier(q, u, q_cos_lm, q_sin_lm, u_cos_lm, u_sin_lm, name,
                          npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1);

        singular_points_classifier(q, u, q_cos_lm, q_sin_lm, u_cos_lm, u_sin_lm, name,
                          npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }
}

void level_area(long double** map, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double current_level = lower_level;

    for (unsigned int i_level = 0; i_level <= nlevel; ++i_level) {

        out_file << std::scientific << current_level << " "
                 << std::scientific << area(map, current_level) << " "
                 << std::endl;

        current_level = current_level + (top_level - lower_level) / nlevel;

    }

    out_file.close();
}

void level_length(long double** map, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double current_level = lower_level;

    for (unsigned int i_level = 0; i_level <= nlevel; ++i_level) {

        out_file << std::scientific << current_level << " "
                 << std::scientific << length(map, current_level) << " "
                 << std::endl;

        current_level = current_level + (top_level - lower_level) / nlevel;

    }

    out_file.close();
}

void level_points(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                  std::string name) {

    long double current_level = lower_level;

    for (unsigned int i_level = 0; i_level <= nlevel; ++i_level) {

        level_points_classifier(map_x, map_y, cos_lm, sin_lm, current_level, name);

        current_level = current_level + (top_level - lower_level) / nlevel;

    }
}

void level_singular_points(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                           long double** u_cos_lm, long double** u_sin_lm, std::string name) {
    long double current_level = lower_level;

    for (unsigned int i_level = 0; i_level <= nlevel; ++i_level) {

        level_singular_points_classifier(q, u, q_cos_lm, q_sin_lm, u_cos_lm, u_sin_lm, current_level, name);

        current_level = current_level + (top_level - lower_level) / nlevel;

    }

}
