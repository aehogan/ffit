#include <axilrod.h>
#include <vector>
#include <cmath>
#include <stdio.h>

double axilrod(System &system, Parameters &parameters)
{
    double energy = 0.;
    vector<Atom> &atoms = system.atoms;
    const int N = atoms.size();
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
        {
            for (int k = j+1; k < N; k++)
            {

                if (atoms[i].mol == atoms[j].mol || atoms[i].mol == atoms[k].mol || atoms[j].mol == atoms[k].mol)
                    continue;

                //printf("---\n");
                
                double dx1, dy1, dz1;
                double dx2, dy2, dz2;
                double dx3, dy3, dz3;
                
                dx1 = atoms[i].x - atoms[j].x;
                dy1 = atoms[i].y - atoms[j].y;
                dz1 = atoms[i].z - atoms[j].z;
                
                dx2 = atoms[i].x - atoms[k].x;
                dy2 = atoms[i].y - atoms[k].y;
                dz2 = atoms[i].z - atoms[k].z;
                
                dx3 = atoms[j].x - atoms[k].x;
                dy3 = atoms[j].y - atoms[k].y;
                dz3 = atoms[j].z - atoms[k].z;

                double3 displacement1, displacement2, displacement3;
                displacement1 = system.pbc.min_image(dx1, dy1, dz1);
                displacement2 = system.pbc.min_image(dx2, dy2, dz2);
                displacement3 = system.pbc.min_image(dx3, dy3, dz3);
                
                dx1 = displacement1.x;
                dy1 = displacement1.y;
                dz1 = displacement1.z;
                
                dx2 = displacement2.x;
                dy2 = displacement2.y;
                dz2 = displacement2.z;
                
                dx3 = displacement3.x;
                dy3 = displacement3.y;
                dz3 = displacement3.z;

                const double r2_1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
                const double r2_2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
                const double r2_3 = dx3*dx3 + dy3*dy3 + dz3*dz3;

                if (r2_1 > system.cutoff * system.cutoff ||
                    r2_2 > system.cutoff * system.cutoff ||
                    r2_3 > system.cutoff * system.cutoff)
                    continue;

                const double r_1 = sqrt(r2_1);
                const double r_2 = sqrt(r2_2);
                const double r_3 = sqrt(r2_3);
                const double c9_1 = 3.0/4.0 * atoms[i].alpha * atoms[i].c6;
                const double c9_2 = 3.0/4.0 * atoms[j].alpha * atoms[j].c6;
                const double c9_3 = 3.0/4.0 * atoms[k].alpha * atoms[k].c6;
                const double c9 = pow(c9_1 * c9_2 * c9_3, 1.0/3.0);
                
                const double cos_1 = ( dx1*dx2 + dy1*dy2 + dz1*dz2 ) / (r_1 * r_2);
                const double cos_2 = ( dx1*dx3 + dy1*dy3 + dz1*dz3 ) / (r_1 * r_3);
                const double cos_3 = ( dx2*dx3 + dy2*dy3 + dz2*dz3 ) / (r_2 * r_3);
                
                //printf("c6 %f %f %f\n", atoms[i].c6, atoms[j].c6, atoms[k].c6);
                //printf("c9 %f\n", c9);
                //printf("r %f %f %f\n", r_1, r_2, r_3);
                //printf("cos %f %f %f\n", cos_1, cos_2, cos_3);
                
                if (parameters.force_fitting_on)
                {
                    // TODO
                }

                if (parameters.energy_fitting_on)
                {
                    energy += c9 * (1 + 3 * cos_1 * cos_2 * cos_3) / (r_1*r_1*r_1 * r_2*r_2*r_2 * r_3*r_3*r_3);
                    //printf("energy %f\n", c9 * (1 + 3 * cos_1 * cos_2 * cos_3) / (r_1*r_1*r_1 * r_2*r_2*r_2 * r_3*r_3*r_3));
                }
            }
        }
    }

    return energy;
}
