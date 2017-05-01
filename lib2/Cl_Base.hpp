#pragma once

class Cl_Base {

private:

    unsigned int n_mod;
    long double* array;

public:

    Cl_Base(unsigned int n_mod);

    virtual ~Cl_Base();

    unsigned int n_mod_value() const;

    long double at(unsigned int l) const;

    void to(unsigned int l, long double cl);

    void fout(std::string name) const;

    void fin(std::string name) const;

};
