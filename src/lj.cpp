#include <lj.h>
#include <vector>
#include <cmath>

double lj(System &system, Parameters &parameters)
{
    double energy = 0.;
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    for(int i = 0; i < N; i++)
    {
        //printf("i %d\n",i);
        for(int j = i+1; j < N; j++)
        {
            //printf("j %d\nN %d\n",j,N);
            //printf("what %d\n",j<N);
            //printf("pair\n---\n");
            //printf("%d %d\n",i,j);
            //printf("%d %d\n",atoms[i].type,atoms[j].type);
            //printf("%d\n",parameters.type_map[string("H")]);
            //printf("%d %d\n",atoms[i].mol,atoms[j].mol);
            //printf("%f %f\n",atoms[i].lj_eps,atoms[j].lj_eps);
            //printf("%f %f\n",atoms[i].lj_sig,atoms[j].lj_sig);

            if (atoms[i].mol == atoms[j].mol)
            {
                //printf("same mol\n\n");
                continue;
            }
            if (atoms[i].lj_eps == 0. || atoms[j].lj_eps == 0. || atoms[i].lj_sig == 0. || atoms[j].lj_sig == 0.)
            {
                //printf("zero lj param(s)\n\n");
                continue;
            }
            
            double dx, dy, dz;
            
            dx = atoms[i].x - atoms[j].x;
            dy = atoms[i].y - atoms[j].y;
            dz = atoms[i].z - atoms[j].z;
            //printf("dx %f %f %f\n",dx,dy,dz);
            double3 displacement;
            displacement = system.pbc.min_image(dx,dy,dz);
            
            dx = displacement.x;
            dy = displacement.y;
            dz = displacement.z;
            //printf("dx %f %f %f\n",dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            //printf("r2 cutoff %f %f\n",r2,system.cutoff);
            if (r2 > system.cutoff * system.cutoff)
                continue;
            
            const double eps = sqrt( atoms[i].lj_eps * atoms[j].lj_eps );
            const double sig = 0.5 * ( atoms[i].lj_sig + atoms[j].lj_sig );
            
            const double r = sqrt(r2);
            
            if (parameters.force_fitting_on)
            {
                const double pre = 24. * eps / r * ( 2. * pow(sig/r,12) - pow(sig/r,6) );
                
                atoms[i].fx += dx*pre;
                atoms[i].fy += dy*pre;
                atoms[i].fz += dz*pre;
                
                atoms[j].fx -= dx*pre;
                atoms[j].fy -= dy*pre;
                atoms[j].fz -= dz*pre;
            }

            if (parameters.energy_fitting_on)
                energy += 4. * eps * ( pow(sig/r,12) - pow(sig/r,6) );
            //printf("r2 %f\nr %f\neps %f\nsig %f\n",r2,r,eps,sig);
            //printf("energy %f\n\n",energy);
        }
    }
    //printf("---\n");
    //printf("lj e %15.6f\n",energy);
    return energy;
}

