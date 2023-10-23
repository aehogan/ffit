#pragma once
#include <vector>
#include <string>
#include <map>

using namespace std;

struct double3 {
    double x, y, z;
};

class Atom {
    public:
        unsigned int type;
        unsigned int mol;
        unsigned int id;
        double x, y, z;
        double fx, fy, fz;
        double fx_fit, fy_fit, fz_fit;
        double charge, lj_eps, lj_sig, alpha, c6, c8, c10;
        double3 mu, ef_static, ef_induced;
};

class PBC {
    public:
        void load_abc(double,double,double,double,double,double);
        double3 min_image(double,double,double) const;
    private:
        double basisx1, basisx2, basisx3;
        double basisy1, basisy2, basisy3;
        double basisz1, basisz2, basisz3;
        double volume, inverse_volume;
        double recip_basisx1, recip_basisx2, recip_basisx3;
        double recip_basisy1, recip_basisy2, recip_basisy3;
        double recip_basisz1, recip_basisz2, recip_basisz3;
};

class Parameters {
    public:
        // parameters
        map<int, double> lj_sig;            // atom type   -> sigma
        map<int, double> lj_eps;            // atom type   -> epsilon
        map<int, double> c6;                // atom type   -> c6
        map<int, double> c8;                // atom type   -> c8
        map<int, double> c10;               // atom type   -> c10
        map<int, double> alpha;             // atom type   -> alpha
        map<int, double> xdm_m1;            // atom type   -> m1
        map<int, double> xdm_m2;            // atom type   -> m2
        map<int, double> xdm_m3;            // atom type   -> m3
        map<string, int> type_map;          // atom string -> atom type
        map<int, string> type_map_reverse;  // atom type   -> atom string
        double wolf_alpha;

        // options
        double cutoff, max_energy;
        unsigned short int axilrod_on, lj_on, phahst_on, es_on, pol_on;
        unsigned short int force_fitting_on, energy_fitting_on, damp_dispersion;
        unsigned short int fit_alpha_on, xdm_on, disp_coeff_extrapolation;
        unsigned int steps, output_freq;
        vector<string> do_not_fit;

        Parameters();
};

class System {
    public:
        // atoms and molecules
        vector<Atom> atoms;
        vector<vector<Atom*>> molecules;
        PBC pbc;

        // internal data
        double cutoff, lj_energy, es_energy, pol_energy, disp_energy, total_energy, fit_energy;

        System();
};

