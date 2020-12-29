/*
A code for computing the non-linear matter bispectrum (BS) based on BiHalofit (arXiv:1911.07886)
written by Ryuichi Takahashi in Dec 2019
a bug fixed in Mar 2020 (Thanks Sven Heydenreich for indicating it)

This is a C++ object-oriented implementation of the original BiHalofit code
such that it can be easily integrated into a c++ project as a class object with easily accessible methods

The default bihalofit class is created with cosmological parameters that are consistent with the Planck 2015 LambdaCDM model.
If you want to initialise the class with another cosmological model, please use the parameterised constructor.

You can also provide a table of linear P(k) (1st & 2nd columns should be k[h/Mpc] & P(k)[(Mpc/h)^3] respectively).
If you give a P(k) table, then the code uses it.
If not, the code uses the Eisenstein & Hu (1999) fitting formula as the linear P(k).

*/

#ifndef BIHALOFIT_HPP
#define BIHALOFIT_HPP

namespace
{
    #define n_data_max 20000
}

class bihalofit
{
public:
    bihalofit(); // default constructor with Planck 2015 cosmology
    bihalofit(double Omega_b, double Omega_cdm, double h, double sigma8, double n_s, double w0); // parameterised constructor for flat wCDM model
    ~bihalofit();

    double omb, omc, h, sigma8, ns, w;
    double om, ow, norm;
    double k_data[n_data_max], pk_data[n_data_max];
    int n_data;

    void compute_Pk_norm(); // P(k) normalisation using sigma8

    // linear pk
    void load_pk_data(char* filename); // load a table of linear P(k) (1st & 2nd columns should be k[h/Mpc] & P(k)[(Mpc/h)^3] respectively).
    double linear_pk_data(double k);
    double linear_pk_eh(double k);
    double linear_pk(double k);

    // linear growth factor
    double lgr_func(int j, double x, double y[2]);
    double lgr(double z);
    double calc_D1(double z);

    // sigma
    double window(double x, int i);
    double sigmam(double r, int j); // Eq.(B1) 2nd part
    double calc_r_sigma(double z); // Eq.(B1) 1st part

    // 3D bispectrum computation formulae
    double F2_tree(double k1, double k2, double k3);
    double bispec_tree(double k1, double k2, double k3, double z);
    double F2(double k1, double k2, double k3, double z);
    double bispec(double k1, double k2, double k3, double z);
    double baryon_ratio(double k1, double k2, double k3, double z);

};

#endif // BIHALOFIT_HPP
