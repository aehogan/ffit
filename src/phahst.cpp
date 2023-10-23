#include <phahst.h>
#include <vector>
#include <cmath>

double phahst(System &system, Parameters &parameters)
{
    double energy = 0.;
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
        {

            if (atoms[i].mol == atoms[j].mol)
                continue;

            if (atoms[i].lj_eps == 0. || atoms[j].lj_eps == 0. || atoms[i].lj_sig == 0. || atoms[j].lj_sig == 0.)
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

            if (r2 > system.cutoff * system.cutoff)
                continue;

            const double r = sqrt(r2);
            const double r4 = r2*r2;
            const double r6 = r4*r2;
            const double r8 = r4*r4;
            const double r10 = r4*r6;
            const double r12 = r6*r6;

            const double c6 = sqrt( atoms[i].c6 * atoms[j].c6 );
            const double c8 = sqrt( atoms[i].c8 * atoms[j].c8 );
            const double c10 = sqrt( atoms[i].c10 * atoms[j].c10 );

            const double rho = 0.5 * ( atoms[i].lj_sig + atoms[j].lj_sig );
            const double beta = 2.0 * atoms[i].lj_eps * atoms[j].lj_eps / ( atoms[i].lj_eps + atoms[j].lj_eps );

            
            if (parameters.force_fitting_on)
            {
                const double pre = 0.;
                
                atoms[i].fx += dx*pre;
                atoms[i].fy += dy*pre;
                atoms[i].fz += dz*pre;
                
                atoms[j].fx -= dx*pre;
                atoms[j].fy -= dy*pre;
                atoms[j].fz -= dz*pre;
            }

            if (parameters.energy_fitting_on)
            {
                double repulsion = 596.725194095 * 1.0 / beta * exp(-beta * (r - rho));

                if (parameters.damp_dispersion)
                    energy += -tt_damping(6, beta * r) * c6 / r6 - tt_damping(8, beta * r) * c8 / r8 - tt_damping(10, beta * r) * c10 / r10 + repulsion;
                else
                    energy += -c6 / r6 - c8 / r8 - c10 / r10 + repulsion;

                //printf("%25.15f %25.15f %25.15f %25.15f %25.15f\n",repulsion,-tt_damping(6, beta * r) * c6 / r6,- tt_damping(8, beta * r) * c8 / r8, - tt_damping(10, beta * r) * c10 / r10, energy);
            }
        }
    }

    return energy;
}

double phahst_dispersion_only(System &system, Parameters &parameters)
{
    double energy = 0.;
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
        {

            if (atoms[i].mol == atoms[j].mol)
                continue;

            if (atoms[i].lj_eps == 0. || atoms[j].lj_eps == 0. || atoms[i].lj_sig == 0. || atoms[j].lj_sig == 0.)
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

            if (r2 > system.cutoff * system.cutoff)
                continue;

            const double r = sqrt(r2);
            const double r4 = r2*r2;
            const double r6 = r4*r2;
            const double r8 = r4*r4;
            const double r10 = r4*r6;
            const double r12 = r6*r6;

            const double c6 = sqrt( atoms[i].c6 * atoms[j].c6 );
            const double c8 = sqrt( atoms[i].c8 * atoms[j].c8 );
            const double c10 = sqrt( atoms[i].c10 * atoms[j].c10 );

            const double beta = 2.0 * atoms[i].lj_eps * atoms[j].lj_eps / ( atoms[i].lj_eps + atoms[j].lj_eps );

            if (parameters.damp_dispersion)
                energy += -tt_damping(6, beta * r) * c6 / r6 - tt_damping(8, beta * r) * c8 / r8 - tt_damping(10, beta * r) * c10 / r10;
            else
                energy += -c6 / r6 - c8 / r8 - c10 / r10;

        }
    }

    return energy;
}

double factorial(int n) {
    int i;
    double fac = 1.0;
    for (i = 2; i <= n; i++)
        fac *= i;
    return fac;
}

double tt_damping(int n, double br) {
    double sum = 1.0, running_br = br;
    int i;
    for (i = 1; i <= n; i++) {
        sum += running_br / factorial(i);
        running_br *= br;
    }

    const double result = 1.0 - exp(-br) * sum;

    if (result > 0.000000001)
        return result;
    else
        return 0.0; /* This is so close to zero lets just call it zero to avoid rounding error and the simulation blowing up */
}

