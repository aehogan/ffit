#include <system.h>
#include <cmath>

Parameters::Parameters()
{
    // default parameters
    wolf_alpha = 0.;

    // default options
    cutoff = 1000.;
    max_energy = 1000.;
    lj_on = 1;
    phahst_on = 0;
    es_on = 1;
    pol_on = 1;
    force_fitting_on = 0;
    energy_fitting_on = 1;
    damp_dispersion = 1;
    steps = 1000000;
    output_freq = 100;
    fit_alpha_on = 0;
    xdm_on = 0;
}

void PBC::load_abc(double a, double b, double c, double alpha, double beta, double gamma)
{   
    basisx1 = a;
    basisx2 = 0.;
    basisx3 = 0.;

    basisy1 = b * cos(M_PI / 180.0 * gamma);
    basisy2 = b * sin(M_PI / 180.0 * gamma);
    basisy3 = 0.;

    basisz1 = c * cos(M_PI / 180.0 * beta);
    basisz2 = ((b * c * cos(M_PI / 180.0 * alpha)) - (basisy1 * basisz1)) / basisy2;
    basisz3 = sqrt(c * c - basisz1 * basisz1 - basisz2 * basisz2);

    volume = basisx1 * (basisy2 * basisz3 - basisy3 * basisz2);
    volume += basisx2 * (basisy3 * basisz1 - basisy1 * basisz3);
    volume += basisx3 * (basisy1 * basisz2 - basisy2 * basisz1);

    inverse_volume = 1.0 / volume;

    recip_basisx1 = inverse_volume * (basisy2 * basisz3 - basisy3 * basisz2);
    recip_basisx2 = inverse_volume * (basisx3 * basisz2 - basisx2 * basisz3);
    recip_basisx3 = inverse_volume * (basisx2 * basisy3 - basisx3 * basisy2);

    recip_basisy1 = inverse_volume * (basisy3 * basisz1 - basisy1 * basisz3);
    recip_basisy2 = inverse_volume * (basisx1 * basisz3 - basisx3 * basisz1);
    recip_basisy3 = inverse_volume * (basisx3 * basisy1 - basisx1 * basisy3);

    recip_basisz1 = inverse_volume * (basisy1 * basisz2 - basisy2 * basisz1);
    recip_basisz2 = inverse_volume * (basisx2 * basisz1 - basisx1 * basisz2);
    recip_basisz3 = inverse_volume * (basisx1 * basisy2 - basisx2 * basisy1);
}

double3 PBC::min_image(double dx, double dy, double dz) const
{
    double imgx = 0., imgy = 0., imgz = 0.;

    imgx += recip_basisx1 * dx;
    imgx += recip_basisy1 * dy;
    imgx += recip_basisz1 * dz;
    imgx = round(imgx);

    imgy += recip_basisx2 * dx;
    imgy += recip_basisy2 * dy;
    imgy += recip_basisz2 * dz;
    imgy = round(imgy);

    imgz += recip_basisx3 * dx;
    imgz += recip_basisy3 * dy;
    imgz += recip_basisz3 * dz;
    imgz = round(imgz);

    double dix = 0., diy = 0., diz = 0.;

    dix += basisx1 * imgx;
    dix += basisy1 * imgy;
    dix += basisz1 * imgz;

    diy += basisx2 * imgx;
    diy += basisy2 * imgy;
    diy += basisz2 * imgz;

    diz += basisx3 * imgx;
    diz += basisy3 * imgy;
    diz += basisz3 * imgz;

    double3 result;
    
    result.x = dx - dix;
    result.y = dy - diy;
    result.z = dz - diz;

    return result;
}

System::System()
{
    // internal data
    lj_energy = 0.;
    es_energy = 0.;
    pol_energy = 0.;
    disp_energy = 0.;
}

