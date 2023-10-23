#include <pol.h>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>

#define OneOverSqrtPi 0.56418958354

void init_dipoles(System &system, const Parameters &parameters)
{
    const double gamma = 1.03;
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    for(int i = 0; i < N; i++)
    {
        atoms[i].mu.x = atoms[i].alpha * atoms[i].ef_static.x * gamma;
        atoms[i].mu.y = atoms[i].alpha * atoms[i].ef_static.y * gamma;
        atoms[i].mu.z = atoms[i].alpha * atoms[i].ef_static.z * gamma;
    }
}

void calc_e_static(System &system, const Parameters &parameters)
{
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    const double a = parameters.wolf_alpha;
    const double cutoffterm = (erfc(a * system.cutoff) * 1./system.cutoff * 1./system.cutoff + 2.0 * a * OneOverSqrtPi * exp(-a * a * system.cutoff * system.cutoff) * 1./system.cutoff);
    
    for (int i = 0; i < N; i++)
    {
        atoms[i].ef_static.x = 0.;
        atoms[i].ef_static.y = 0.;
        atoms[i].ef_static.z = 0.;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            if (atoms[i].mol == atoms[j].mol)
                continue;
            
            double dx, dy, dz;
            
            dx = atoms[i].x - atoms[j].x;
            dy = atoms[i].y - atoms[j].y;
            dz = atoms[i].z - atoms[j].z;
            
            double3 displacement;
            displacement = system.pbc.min_image(dx,dy,dz);
            
            dx = displacement.x;
            dy = displacement.y;
            dz = displacement.z;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            if (r2 > system.cutoff*system.cutoff)
                continue;
            
            const double r = sqrt(r2);
            const double ir = 1. / r;

            if (a == 0.) {

                atoms[i].ef_static.x += atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dx * ir;
                atoms[i].ef_static.y += atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dy * ir;
                atoms[i].ef_static.z += atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dz * ir;
                atoms[j].ef_static.x -= atoms[i].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dx * ir;
                atoms[j].ef_static.y -= atoms[i].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dy * ir;
                atoms[j].ef_static.z -= atoms[i].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dz * ir;

            } else {

                const double bigmess = (erfc(a * r) * ir * ir + 2.0 * a * OneOverSqrtPi * exp(-a * a * r * r) * ir);
                atoms[i].ef_static.x += atoms[j].charge * (bigmess - cutoffterm) * dx * ir;
                atoms[i].ef_static.y += atoms[j].charge * (bigmess - cutoffterm) * dy * ir;
                atoms[i].ef_static.z += atoms[j].charge * (bigmess - cutoffterm) * dz * ir;
                atoms[j].ef_static.x -= atoms[i].charge * (bigmess - cutoffterm) * dx * ir;
                atoms[j].ef_static.y -= atoms[i].charge * (bigmess - cutoffterm) * dy * ir;
                atoms[j].ef_static.z -= atoms[i].charge * (bigmess - cutoffterm) * dz * ir;

            }
        }
    }
}

bool pair_cmp(const pair<int, double>& a, const pair<int, double>& b)
{
    return a.second > b.second;
}

void make_ranked_array(const System &system, const Parameters &parameters, vector<int> &ranked_array)
{
    const vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    vector<pair<int, double>> pairs;

    const double rmin = 4.0;

    // calculate rank metrics
    for (int i = 0; i < N; i++)
    {
        double rank_metric = 0.;
        if (atoms[i].alpha != 0.)
        {
            for (int j = 0; j < N; j++)
            {
                if (i == j)
                    continue;

                double dx, dy, dz;
            
                dx = atoms[i].x - atoms[j].x;
                dy = atoms[i].y - atoms[j].y;
                dz = atoms[i].z - atoms[j].z;
                
                double3 displacement;
                displacement = system.pbc.min_image(dx,dy,dz);
                
                dx = displacement.x;
                dy = displacement.y;
                dz = displacement.z;
                
                const double r2 = dx*dx + dy*dy + dz*dz;

                if (r2 > rmin*rmin)
                    continue;

                rank_metric += atoms[i].alpha;
            }
        }
        pair<int, double> pair;
        pair.first = i;
        pair.second = rank_metric;
        pairs.push_back(pair);
    }

    // sort
    sort(pairs.begin(), pairs.end(), pair_cmp);

    for (auto &pair : pairs)
        ranked_array.push_back(pair.first);

}

void iterate_pol(System &system, const Parameters &parameters, const vector<int> &ranked_array)
{
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();

    for (int ii = 0; ii < N; ii++)
    {
        const int i = ranked_array[ii];

        atoms[i].ef_induced.x = 0.;
        atoms[i].ef_induced.y = 0.;
        atoms[i].ef_induced.z = 0.;

        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;

            double dx, dy, dz;
            
            dx = atoms[i].x - atoms[j].x;
            dy = atoms[i].y - atoms[j].y;
            dz = atoms[i].z - atoms[j].z;
                
            double3 displacement;
            displacement = system.pbc.min_image(dx,dy,dz);
              
            dx = displacement.x;
            dy = displacement.y;
            dz = displacement.z;
                
            const double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > system.cutoff*system.cutoff)
                continue;

            const double l = 2.1304;

            const double r = sqrt(r2);
            const double ir = 1./r;
            const double ir2 = 1./r2;
            const double ir3 = ir*ir2;
            const double ir5 = ir2*ir3;

            const double explr = exp(-l * r);
            const double damp1 = 1.0 - explr * (0.5 * l * l * r2 + l * r + 1.0);
            const double damp2 = damp1 - explr * (l * l * l * r2 * r / 6.0);

            atoms[i].ef_induced.x -= -3.0 * dx * dx * damp2 * ir5 * atoms[j].mu.x;
            atoms[i].ef_induced.x -= -3.0 * dx * dy * damp2 * ir5 * atoms[j].mu.y;
            atoms[i].ef_induced.x -= -3.0 * dx * dz * damp2 * ir5 * atoms[j].mu.z;

            atoms[i].ef_induced.y -= -3.0 * dy * dx * damp2 * ir5 * atoms[j].mu.x;
            atoms[i].ef_induced.y -= -3.0 * dy * dy * damp2 * ir5 * atoms[j].mu.y;
            atoms[i].ef_induced.y -= -3.0 * dy * dz * damp2 * ir5 * atoms[j].mu.z;

            atoms[i].ef_induced.z -= -3.0 * dz * dx * damp2 * ir5 * atoms[j].mu.x;
            atoms[i].ef_induced.z -= -3.0 * dz * dy * damp2 * ir5 * atoms[j].mu.y;
            atoms[i].ef_induced.z -= -3.0 * dz * dz * damp2 * ir5 * atoms[j].mu.z;

            atoms[i].ef_induced.x -= damp1 * ir3 * atoms[j].mu.x;
            atoms[i].ef_induced.y -= damp1 * ir3 * atoms[j].mu.y;
            atoms[i].ef_induced.z -= damp1 * ir3 * atoms[j].mu.z;

        }

        atoms[i].mu.x = atoms[i].alpha * ( atoms[i].ef_static.x + atoms[i].ef_induced.x );
        atoms[i].mu.y = atoms[i].alpha * ( atoms[i].ef_static.y + atoms[i].ef_induced.y );
        atoms[i].mu.z = atoms[i].alpha * ( atoms[i].ef_static.z + atoms[i].ef_induced.z );
    }
}

double pol(System &system, Parameters &parameters)
{
    calc_e_static(system,parameters);
    init_dipoles(system,parameters);

    vector<int> ranked_array;
    make_ranked_array(system,parameters,ranked_array);

    const int iterations = 3;

    for (int i = 0; i < iterations; i++)
      iterate_pol(system,parameters,ranked_array);

    double energy = 0.;

    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();

    for (int i = 0; i < N; i++)
    {
        energy += -0.5 * atoms[i].ef_static.x * atoms[i].mu.x;
        energy += -0.5 * atoms[i].ef_static.y * atoms[i].mu.y;
        energy += -0.5 * atoms[i].ef_static.z * atoms[i].mu.z;
    }

/*
    printf("mu MOLECULE ATOM * DIPOLES * STATIC * INDUCED * pot/atom -0.5*mu*E_s\n");
    for (int i = 0; i < N; i++)
        printf("mu %4d %4d * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %lf %lf\n",
            atoms[i].mol, atoms[i].id,
            atoms[i].mu.x, atoms[i].mu.y, atoms[i].mu.z,
            atoms[i].ef_static.x, atoms[i].ef_static.y, atoms[i].ef_static.z,
            atoms[i].ef_induced.x, atoms[i].ef_induced.y, atoms[i].ef_induced.z,
            energy / N, -0.5 * atoms[i].mu.x * atoms[i].ef_static.x);
*/

    return energy;
}

