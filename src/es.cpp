#include <es.h>
#include <vector>
#include <cmath>

#define OneOverSqrtPi 0.56418958354

double es(System &system, Parameters &parameters)
{
    double energy = 0.;
    vector<Atom> &atoms = system.atoms;
    int N = atoms.size();
    const double a = parameters.wolf_alpha;
    const double cutoffterm = (erfc(a * system.cutoff) * 1./system.cutoff * 1./system.cutoff + 2.0 * a * OneOverSqrtPi * exp(-a * a * system.cutoff * system.cutoff) * 1./system.cutoff);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
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
                if (parameters.force_fitting_on)
                {
                    atoms[i].fx += atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dx * ir;
                    atoms[i].fy += atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dy * ir;
                    atoms[i].fz += atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dz * ir;
                    atoms[j].fx -= atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dx * ir;
                    atoms[j].fy -= atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dy * ir;
                    atoms[j].fz -= atoms[i].charge * atoms[j].charge * (ir * ir - 1./system.cutoff * 1./system.cutoff) * dz * ir;
                }
                if (parameters.energy_fitting_on)
                    energy += atoms[i].charge * atoms[j].charge * ( ir - 1./system.cutoff);
            } else {
                if (parameters.force_fitting_on)
                {
                    const double bigmess = (erfc(a * r) * ir * ir + 2.0 * a * OneOverSqrtPi * exp(-a * a * r * r) * ir);
                    atoms[i].fx += atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dx * ir;
                    atoms[i].fy += atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dy * ir;
                    atoms[i].fz += atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dz * ir;
                    atoms[j].fx -= atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dx * ir;
                    atoms[j].fy -= atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dy * ir;
                    atoms[j].fz -= atoms[i].charge * atoms[j].charge * (bigmess - cutoffterm) * dz * ir;
                }
                if (parameters.energy_fitting_on)
                    if (a == 0.)
                        energy += atoms[i].charge * atoms[j].charge * ( ir - 1. / system.cutoff);
                    else
                        energy += atoms[i].charge * atoms[j].charge * ( erfc(a * r) * ir - erfc(a * system.cutoff) / system.cutoff);
            }

        }
    }

    if ( a != 0. && parameters.energy_fitting_on )
    {
        const double self_cutoff_term = erfc(a * system.cutoff) / (2. * system.cutoff ) + a * OneOverSqrtPi;
        for (int i = 0; i < N; i++)
        {
            energy += self_cutoff_term * atoms[i].charge * atoms[i].charge;
        }
    }

    return energy;
}

