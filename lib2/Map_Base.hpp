#pragma once

class Map_Base {

private:

    unsigned int n_pix;
    long double** array;

public:

    Map_Base(unsigned int n_pix_input);

    virtual ~Map_Base();

    unsigned int n_pix_value() const;

    long double at(unsigned int i, unsigned int j) const;

    void to(unsigned int i, unsigned int j, long double pix);

    long double at_border(unsigned int i, unsigned int j) const;

    void to_border(unsigned int i, unsigned int j, long double pix);

    void fout(std::string name) const;

    void fin(std::string name) const;
};
