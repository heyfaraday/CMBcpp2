#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

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
const unsigned int monte_n = static_cast<unsigned int>(read_from_file("monte_n"));

const long double long_npix = static_cast<long double>(npix);
const long double long_nmod = static_cast<long double>(nmod);
const long double long_monte_n = static_cast<long double>(monte_n);