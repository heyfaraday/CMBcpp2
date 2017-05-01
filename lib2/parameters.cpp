#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "constants.hpp"

#include "parameters.hpp"

bool is_number(std::string& s);

double read_from_file(const std::string param);

bool is_number(std::string& s) {

    if (s.empty())
        return false;

    bool trigger = false;

    if (s[0] == '.')
        return false;

    for (int i = 0; i < int(s.size()); i ++) {

        char c = s[static_cast<unsigned long>(i)];

        if (c == '.' && !trigger) {
            trigger = !trigger;
            continue;
        } else if (c == '.') {
            return false;
        }

        if (i == 0 && c == '-')
            continue;

        if (!std::isdigit(c))
            return false;
    }

    return true;
}

double read_from_file(const std::string param) {

    std::ifstream infile("input.parameters");
    std::string str;

    if (infile.is_open()) {

        while (getline(infile, str)) {

            std::istringstream iss(str);
            std::vector<std::string> words{std::istream_iterator<std::string>{iss},
                                           std::istream_iterator<std::string>{}};

            if (words.size() < 1)
                continue;

            if (words[0] == param) {

                if (is_number(words[1])) {
                    return atof(words[1].c_str());
                } else {
                    std::cout << "Error converting string to number for parameter " << param
                              << " in the input file.\n";
                    exit (EXIT_FAILURE);
                }
            }
        }
    } else {
        std::cout << "Cannot open the file.\n";
        exit (EXIT_FAILURE);
    }

    std::cout << "Cannot find the given parameter " << param << " in the input file.\n";
    exit (EXIT_FAILURE);
}

const unsigned int npix = static_cast<unsigned int>(read_from_file("npix"));
const unsigned int nmod = static_cast<unsigned int>(read_from_file("nmod"));
const unsigned int nback = static_cast<unsigned int>(read_from_file("nback"));
const unsigned int nmonte = static_cast<unsigned int>(read_from_file("nmonte"));
const unsigned int nring = static_cast<unsigned int>(read_from_file("nring"));
const unsigned int hring = static_cast<unsigned int>(read_from_file("hring"));

const long double long_npix = static_cast<long double>(npix);
const long double long_nmod = static_cast<long double>(nmod);
const long double long_nback = static_cast<long double>(nback);
const long double long_nmonte = static_cast<long double>(nmonte);
const long double long_nring = static_cast<long double>(nring);
const long double long_hring = static_cast<long double>(hring);

const long double map_parameter = 2.0L * PI / long_npix;

const long double lower_level = static_cast<long double>(read_from_file("lower_level"));
const long double top_level = static_cast<long double>(read_from_file("top_level"));
const long double nlevel = static_cast<long double>(read_from_file("nlevel"));
