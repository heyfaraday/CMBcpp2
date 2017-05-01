#pragma once

class Alm_Base {

private:

    unsigned int n_mod;
    long double** array_cos;
    long double** array_sin;

public:

    Alm_Base(unsigned int n_mod);

    virtual ~Alm_Base();

    unsigned int n_mod_value() const;

    long double at_cos(unsigned int l, unsigned int m) const;

    void to_cos(unsigned int l, unsigned int m, long double alm);

    long double at_sin(unsigned int l, unsigned int m) const;

    void to_sin(unsigned int l, unsigned int m, long double alm);

    void fout(std::string name_cos, std::string name_sin) const;

    void fin(std::string name_cos, std::string name_sin) const;

};
