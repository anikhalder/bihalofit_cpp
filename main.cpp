/*
A code for computing the non-linear matter bispectrum (BS) based on BiHalofit (arXiv:1911.07886)
written by Ryuichi Takahashi in Dec 2019
a bug fixed in Mar 2020 (Thanks Sven Heydenreich for indicating it)

This is a C++ object-oriented implementation of the original BiHalofit code
such that it can be easily integrated into a c++ project as a class object with easily accessible methods

The default bihalofit class is created with cosmological parameters that are consistent with the Planck 2015 LambdaCDM model.
If you want to initialise the class with another cosmological model, please use the parameterised constructor.

You can also provide a table of linear linear P(k) (1st & 2nd columns should be k[h/Mpc] & P(k)[(Mpc/h)^3] respectively).
If you give a P(k) table, then the code uses it.
If not, the code uses the Eisenstein & Hu (1999) fitting formula as the linear P(k).

*/


#include <iostream>
#include <bihalofit.hpp>

int main()
{
    // create bihalofit object
    bihalofit bh_obj = bihalofit(); //using default constructor
    //bihalofit bh_obj = bihalofit(0.047, 0.286, 0.7, 0.82, 0.96, -1); //using parameterised constructor

    // can load desired linear power spectrum P(k) at z=0; if not code will use Eisenstein & Hu (1999) fitting formula for P(k)
    //bh_obj.load_pk_data("linear_pk_planck2015.txt");

    double z=0.4;  // redshift
    double k1=1., k2=1.5, k3=2.;   // [h/Mpc]

    // calculating D1, r_sigma & n_eff for a given redshift
    printf("%lf %lf %lf \n", bh_obj.calc_D1(z), bh_obj.calc_r_sigma(z), bh_obj.calc_n_eff(z));

    // non-linear BS w/o baryons, tree-level BS [(Mpc/h)^6] & baryon ratio for given k1,k2,k3 and z
    printf("%lf %lf %lf %lf %lf %lf \n", k1, k2, k3, bh_obj.bispec(k1,k2,k3,z), bh_obj.bispec_tree(k1,k2,k3,z), bh_obj.baryon_ratio(k1,k2,k3,z));

    return 0;
}
