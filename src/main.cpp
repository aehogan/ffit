#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <float.h>

#include <main.h>
#include <system.h>
#include <es.h>
#include <lj.h>
#include <pol.h>
#include <phahst.h>
#include <axilrod.h>

#define DEBUG 0

#define E2REDUCED 408.7816          /* convert from e to sqrt(K*A) */

using namespace std;

int load_systems(vector<System> &systems, Parameters &parameters, vector<string> &data_files)
{
    // check closer for errors here and complain
    unsigned int i = 0;
    for (const auto& data_file : data_files)
    {
        fstream file;
        file.open(data_file,ios::in);
        if (file.is_open())
        {
            string line;
            char cline[1000];
            while (getline(file,line))
            {
                System system;
                system.cutoff = parameters.cutoff;
                strcpy(cline,line.c_str());
                const int N = atoi(strtok(cline," \t"));
                getline(file,line);
                strcpy(cline,line.c_str());
                double energy = atof(strtok(cline," \t"));
                system.fit_energy = energy;
                double a, b, c, alpha, beta, gamma;
                char *token = strtok(NULL," \t");
                if ( token == NULL )
                {
                    a = 1000.;
                    b = 1000.;
                    c = 1000.;
                    alpha = 90.;
                    beta = 90.;
                    gamma = 90.;
                }
                else
                {
                    a = atof(token);
                    b = atof(strtok(NULL," \t"));
                    c = atof(strtok(NULL," \t"));
                    alpha = atof(strtok(NULL," \t"));
                    beta = atof(strtok(NULL," \t"));
                    gamma = atof(strtok(NULL," \t"));
                }
                system.pbc.load_abc(a,b,c,alpha,beta,gamma);
                system.cutoff = min(a,min(b,c))/2.;
                for (int i = 0; i < N; i++)
                {
                    Atom atom;
                    getline(file,line);
                    strcpy(cline,line.c_str());
                    //printf("%s\n",cline);
                    char *token = strtok(cline," \t");
                    atom.type = parameters.type_map[token];
                    if (atom.type == 0)
                    {
                        printf("Atom type %s in %s has no corresponding type_map command in input file\n",token,data_file.c_str());
                        return 1;
                    }
                    atom.mol = atoi(strtok(NULL," \t"));
                    atom.id = i;
                    atom.x = atof(strtok(NULL," \t"));
                    atom.y = atof(strtok(NULL," \t"));
                    atom.z = atof(strtok(NULL," \t"));
                    atom.charge = atof(strtok(NULL," \t"));
                    atom.charge *= E2REDUCED;
                    if (parameters.force_fitting_on)
                    {
                        atom.fx_fit = atof(strtok(NULL," \t"));
                        atom.fy_fit = atof(strtok(NULL," \t"));
                        atom.fz_fit = atof(strtok(NULL," \t"));
                    }
                    else
                    {
                        atom.fx_fit = 0.;
                        atom.fy_fit = 0.;
                        atom.fz_fit = 0.;
                    }
                    atom.alpha = 0.;
                    atom.lj_eps = 0.;
                    atom.lj_sig = 0.;
                    atom.c6 = 0.;
                    atom.c8 = 0.;
                    atom.c10 = 0.;
                    atom.mu.x = 0.;
                    atom.mu.y = 0.;
                    atom.mu.z = 0.;
                    atom.ef_static.x = 0.;
                    atom.ef_static.y = 0.;
                    atom.ef_static.z = 0.;
                    atom.ef_induced.x = 0.;
                    atom.ef_induced.y = 0.;
                    atom.ef_induced.z = 0.;
                    system.atoms.push_back(atom);
                }
                systems.push_back(system);
            }
        }
        else
        {
            printf("Error opening data file %s\n",data_file.c_str());
            return 1;
        }
    }
    return 0;
}

void apply_parameters(vector<System> &systems, Parameters &parameters)
{
    if (parameters.xdm_on)
    {
        for (auto &pair : parameters.type_map)
        {
            int type = pair.second;

            bool do_not_fit_this = false;
            for (auto &name : parameters.do_not_fit)
                if (pair.first == name)
                    do_not_fit_this = true;

            if (do_not_fit_this)
                continue;

            parameters.c6[type] = 0.5 * parameters.alpha[type] * parameters.xdm_m1[type];
            parameters.c6[type] *= 1.8897259886 * 1.8897259886 * 1.8897259886;
            parameters.c6[type] *= 0.021958709 / (3.166811429 * 0.000001);

            parameters.c8[type] = 3. / 2. * parameters.alpha[type] * parameters.xdm_m2[type];
            parameters.c8[type] *= 1.8897259886 * 1.8897259886 * 1.8897259886;
            parameters.c8[type] *= 0.0061490647 / (3.166811429 * 0.000001);

            parameters.c10[type] = 2. * parameters.alpha[type] * parameters.xdm_m3[type] + 21./10. * parameters.alpha[type] * parameters.xdm_m2[type] * parameters.xdm_m2[type] / parameters.xdm_m1[type];
            parameters.c10[type] *= 1.8897259886 * 1.8897259886 * 1.8897259886;
            parameters.c10[type] *= 0.0017219135 / (3.166811429 * 0.000001);
        }
    }

    if (parameters.disp_coeff_extrapolation)
    {
        for (auto &pair : parameters.type_map)
        {
            int type = pair.second;

            bool do_not_fit_this = false;
            for (auto &name : parameters.do_not_fit)
                if (pair.first == name)
                    do_not_fit_this = true;

            if (do_not_fit_this)
                continue;

            parameters.c10[type] = 49./40. * parameters.c8[type] * parameters.c8[type] / parameters.c6[type];
        }
    }

    for (auto &system : systems)
    {
        for (auto &atom : system.atoms)
        {
            atom.lj_eps = parameters.lj_eps[atom.type];
            atom.lj_sig = parameters.lj_sig[atom.type];
            atom.alpha = parameters.alpha[atom.type];
            atom.c6 = parameters.c6[atom.type];
            atom.c8 = parameters.c8[atom.type];
            atom.c10 = parameters.c10[atom.type];
        }
    }
}

double calc_error(vector<System> &systems, Parameters &parameters)
{
    double error = 0.;
    int i = 0;

    if (parameters.force_fitting_on)
    {
        for (const auto &system : systems)
        {
            for (const auto &atom : system.atoms)
            {
                error += pow(atom.fx-atom.fx_fit,2);
                error += pow(atom.fy-atom.fy_fit,2);
                error += pow(atom.fz-atom.fz_fit,2);
                i++;
            }
        }
    }

    if (parameters.energy_fitting_on)
    {
        for (const auto &system : systems)
        {
            double fit_energy, calc_energy;

            fit_energy = system.fit_energy <= 0. ? system.fit_energy : parameters.max_energy * atan(system.fit_energy/parameters.max_energy);
            calc_energy = system.total_energy <= 0. ? system.total_energy : parameters.max_energy * atan(system.total_energy/parameters.max_energy);

            error += pow(calc_energy-fit_energy,2);
            i++;
        }
    }

    if (parameters.lj_on || parameters.phahst_on)
    {
        for (auto &param : parameters.lj_eps)
            if (param.second < 0.)
                return DBL_MAX;
        
        for (auto &param : parameters.lj_sig)
            if (param.second < 0.)
                return DBL_MAX;

        for (auto &param : parameters.c6)
            if (param.second < 0.)
                return DBL_MAX;

        for (auto &param : parameters.c8)
            if (param.second < 0.)
                return DBL_MAX;

        for (auto &param : parameters.c10)
            if (param.second < 0.)
                return DBL_MAX;

        for (auto &param : parameters.alpha)
            if (param.second < 0.)
                return DBL_MAX;
    }

    return error/i;
}

void clear_forces(System &system)
{
    for(auto &atom : system.atoms)
    {
        atom.fx = 0.;
        atom.fy = 0.;
        atom.fz = 0.;
    }
}

void calc_forces(vector<System> &systems, Parameters &parameters)
{
    for (auto &system : systems)
    {
        double energy = 0.;
        clear_forces(system);
        if (parameters.lj_on)
            energy += lj(system, parameters);
        if (parameters.phahst_on)
            energy += phahst(system, parameters);
        if (parameters.es_on)
            energy += es(system, parameters);
        if (parameters.axilrod_on)
            energy += axilrod(system, parameters);
        if (parameters.pol_on && energy <= 10000.)
            energy += pol(system, parameters);
        system.total_energy = energy;
    }
}

double train(vector<System> &systems, Parameters &parameters, bool reset)
{
    const double eps_dx = 0.00001;
    const double eps_alpha = 0.1;
    const double sig_dx = 0.000001;
    const double sig_alpha = 0.5;
    const double alpha_dx = 0.000001;
    const double alpha_alpha = 1.0;
    const double tau = 0.5, c = 0.5;

    static vector<double> adaptive_step_size;
    int param_num = 0;
    
    if (reset)
        adaptive_step_size.clear();

    apply_parameters(systems,parameters);
    calc_forces(systems,parameters);
    double old_error = calc_error(systems,parameters);

    if (parameters.lj_on || parameters.phahst_on)
    {
        for (auto &param : parameters.lj_eps)
        {
            bool do_not_fit_this = false;
            for (auto &name : parameters.do_not_fit)
                if (parameters.type_map_reverse[param.first] == name)
                    do_not_fit_this = true;

            if (do_not_fit_this)
                continue;

            double start_val = param.second;

            param.second += eps_dx;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            double new_error = calc_error(systems,parameters);

            const double slope = ( new_error - old_error ) / eps_dx;
            const double t = c * abs(slope);

            if (param_num >= adaptive_step_size.size())
                adaptive_step_size.push_back(eps_alpha);

            double alpha = slope > 0 ? -adaptive_step_size[param_num] : adaptive_step_size[param_num];
#if DEBUG
            printf("eps start %20.15f %20.15f %20.15f %20.15f %20.15f\n",start_val,error,new_error,slope,alpha);
#endif

            param.second = start_val + alpha;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            new_error = calc_error(systems,parameters);
#if DEBUG
            printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            int linesearch = 0;
            while (old_error - new_error < abs(alpha) * t && abs(alpha) > 1e-1*eps_dx)
            {
                alpha *= tau;

                param.second = start_val + alpha;
                apply_parameters(systems,parameters);
                calc_forces(systems,parameters);
                new_error = calc_error(systems,parameters);
                linesearch++;
#if DEBUG
                printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            }
            old_error = new_error;
            adaptive_step_size[param_num] *= linesearch > 1 ? pow(tau,linesearch-1) : ( linesearch = 0 ? 1./tau : 1. );
            param_num++;
        }

        for (auto &param : parameters.lj_sig)
        {
            bool do_not_fit_this = false;
            for (auto &name : parameters.do_not_fit)
                if (parameters.type_map_reverse[param.first] == name)
                    do_not_fit_this = true;

            if (do_not_fit_this)
                continue;

            double start_val = param.second;

            param.second += sig_dx;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            double new_error = calc_error(systems,parameters);

            const double slope = ( new_error - old_error ) / sig_dx;
            const double t = c * abs(slope);

            if (param_num >= adaptive_step_size.size())
                adaptive_step_size.push_back(sig_alpha);

            double alpha = slope > 0 ? -adaptive_step_size[param_num] : adaptive_step_size[param_num];
#if DEBUG
            printf("sig start %20.15f %20.15f %20.15f %20.15f %20.15f\n",start_val,error,new_error,slope,alpha);
#endif

            param.second = start_val + alpha;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            new_error = calc_error(systems,parameters);
#if DEBUG
            printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            int linesearch = 0;
            while (old_error - new_error < abs(alpha) * t && abs(alpha) > 1e-1*sig_dx)
            {
                alpha *= tau;

                param.second = start_val + alpha;
                apply_parameters(systems,parameters);
                calc_forces(systems,parameters);
                new_error = calc_error(systems,parameters);
                linesearch++;
#if DEBUG
                printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            }
            old_error = new_error;
            adaptive_step_size[param_num] *= linesearch > 1 ? pow(tau,linesearch-1) : ( linesearch = 0 ? 1./tau : 1. );
            param_num++;
        }
    }

    if (parameters.fit_alpha_on)
    {
        for (auto &param : parameters.alpha)
        {
            bool do_not_fit_this = false;
            for (auto &name : parameters.do_not_fit)
                if (parameters.type_map_reverse[param.first] == name)
                    do_not_fit_this = true;

            if (do_not_fit_this)
                continue;

            double start_val = param.second;

            param.second += alpha_dx;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            double new_error = calc_error(systems,parameters);

            const double slope = ( new_error - old_error ) / alpha_dx;
            const double t = c * abs(slope);

            if (param_num >= adaptive_step_size.size())
                adaptive_step_size.push_back(alpha_alpha);

            double alpha = slope > 0 ? -adaptive_step_size[param_num] : adaptive_step_size[param_num];
#if DEBUG
            printf("alpha start %20.15f %20.15f %20.15f %20.15f %20.15f\n",start_val,error,new_error,slope,alpha);
#endif

            param.second = start_val + alpha;
            apply_parameters(systems,parameters);
            calc_forces(systems,parameters);
            new_error = calc_error(systems,parameters);
#if DEBUG
            printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            int linesearch = 0;
            while (old_error - new_error < abs(alpha) * t && abs(alpha) > 1e-1*alpha_dx)
            {
                alpha *= tau;

                param.second = start_val + alpha;
                apply_parameters(systems,parameters);
                calc_forces(systems,parameters);
                new_error = calc_error(systems,parameters);
                linesearch++;
#if DEBUG
                printf("%20.15f %20.15f %20.15f %20.15f %20.15f\n",alpha,error,new_error,error - new_error,abs(alpha) * t);
#endif
            }
            old_error = new_error;
            adaptive_step_size[param_num] *= linesearch > 1 ? pow(tau,linesearch-1) : ( linesearch = 0 ? 1./tau : 1. );
            param_num++;
        }
    }

    return old_error;
}

int load_options(vector<System> &systems, Parameters &parameters, string &input_file)
{
    vector<string> data_files;

    fstream file;
    file.open(input_file,ios::in);
    if (file.is_open())
    {
        string line;
        char cline[1000];
        while (getline(file,line))
        {
            if (line[0] == '\n' || line[0] == '#' || line[0] == '!' || line[0] == '\0')
                continue;

            strcpy(cline,line.c_str());
            char *token = strtok(cline," \t");
            if (strcmp(token,"data_file") == 0)
            {
                token = strtok(NULL," \t");
                string data_file_name = string(token);
                data_files.push_back(data_file_name);
            }
            else if (strcmp(token,"steps") == 0)
            {
                token = strtok(NULL," \t");
                parameters.steps = atoi(token);
            }
            else if (strcmp(token,"output_freq") == 0)
            {
                token = strtok(NULL," \t");
                parameters.output_freq = atoi(token);
            }
            else if (strcmp(token,"cutoff") == 0)
            {
                token = strtok(NULL," \t");
                parameters.cutoff = atof(token);
            }
            else if (strcmp(token,"max_energy") == 0)
            {
                token = strtok(NULL," \t");
                parameters.max_energy = atof(token);
            }
            else if (strcmp(token,"lj") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.lj_on = 1;
                    parameters.phahst_on = 0;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.lj_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"phahst") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.lj_on = 0;
                    parameters.phahst_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.phahst_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"es") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.es_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.es_on = 0;
                }
                else
                    return 1;
            }
            else if (strcmp(token,"axilrod") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.axilrod_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.axilrod_on = 0;
                }
                else
                    return 1;
            }
            else if (strcmp(token,"pol") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.pol_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.pol_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"force_fitting") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.force_fitting_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.force_fitting_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"energy_fitting") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.energy_fitting_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.energy_fitting_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"type_map") == 0)
            {
                token = strtok(NULL," \t");
                string name = string(token);
                token = strtok(NULL," \t");
                int type = atoi(token);
                if (type == 0)
                {
                    printf("Can not use 0 in type_map\n%s\n",line.c_str());
                    return 1;
                }
                parameters.type_map.insert(pair<string,int>(name,type));
            }
            else if (strcmp(token,"lj_sig") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double sig = atof(token);
                parameters.lj_sig.insert(pair<int,double>(type,sig));
            }
            else if (strcmp(token,"lj_eps") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double eps = atof(token);
                parameters.lj_eps.insert(pair<int,double>(type,eps));
            }
            else if (strcmp(token,"c6") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double c6 = atof(token);
                c6 *= 0.021958709 / (3.166811429 * 0.000001);
                parameters.c6.insert(pair<int,double>(type,c6));
            }
            else if (strcmp(token,"c8") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double c8 = atof(token);
                c8 *= 0.0061490647 / (3.166811429 * 0.000001);
                parameters.c8.insert(pair<int,double>(type,c8));
            }
            else if (strcmp(token,"c10") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double c10 = atof(token);
                c10 *= 0.0017219135 / (3.166811429 * 0.000001);
                parameters.c10.insert(pair<int,double>(type,c10));
            }
            else if (strcmp(token,"alpha") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double alpha = atof(token);
                parameters.alpha.insert(pair<int,double>(type,alpha));
            }
            else if (strcmp(token,"fit_alpha") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.fit_alpha_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.fit_alpha_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"xdm") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.xdm_on = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.xdm_on = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"disp_coeff_extrapolation") == 0)
            {
                token = strtok(NULL," \t");
                if (strcmp(token,"on") == 0)
                {
                    parameters.disp_coeff_extrapolation = 1;
                }
                else if (strcmp(token,"off") == 0)
                {
                    parameters.disp_coeff_extrapolation = 0;
                }
                else
                {
                    printf("Error on line:\n%s\n",line.c_str());
                    return 1;
                }
            }
            else if (strcmp(token,"xdm_m1") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double xdm = atof(token);
                parameters.xdm_m1.insert(pair<int,double>(type,xdm));
            }
            else if (strcmp(token,"xdm_m2") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double xdm = atof(token);
                parameters.xdm_m2.insert(pair<int,double>(type,xdm));
            }
            else if (strcmp(token,"xdm_m3") == 0)
            {
                token = strtok(NULL," \t");
                int type = atoi(token);
                token = strtok(NULL," \t");
                double xdm = atof(token);
                parameters.xdm_m3.insert(pair<int,double>(type,xdm));
            }
            else if (strcmp(token,"do_not_fit") == 0)
            {
                token = strtok(NULL," \t");
                while (token != NULL)
                {
                    string name = string(token);
                    parameters.do_not_fit.push_back(name);
                    token = strtok(NULL," \t");
                }
            }
            else
            {
                printf("Error on line:\n%s\n",line.c_str());
                return 1;
            }

        }
    }
    else
    {
        printf("Error opening input file\n");
        return 1;
    }

    // TODO phahst and lj mutually exclusive

    printf("\nSelected options\n\n");
    printf("steps: %d\n",parameters.steps);
    printf("output_freq: %d\n",parameters.output_freq);
    printf("max_energy: %f\n",parameters.max_energy);
    if (parameters.lj_on)
        printf("lj: on\n");
    else
        printf("lj: off\n");
    if (parameters.axilrod_on)
        printf("axilrod: on\n");
    else
        printf("axilrod: off\n");
    if (parameters.phahst_on)
        printf("phahst: on\n");
    else
        printf("phahst: off\n");
    if (parameters.es_on)
        printf("es: on\n");
    else
        printf("es: off\n");
    if (parameters.pol_on)
        printf("pol: on\n");
    else
        printf("pol: off\n");
    if (parameters.force_fitting_on)
        printf("force fitting: on\n");
    else
        printf("force fitting: off\n");
    if (parameters.energy_fitting_on)
        printf("energy fitting: on\n");
    else
        printf("energy fitting: off\n");

    printf("\nLoading systems ...\n");
    if (load_systems(systems,parameters,data_files)!=0)
    {
        printf("Error in loading data files\n");
        return 1;
    }

    for (auto i = parameters.type_map.begin(); i != parameters.type_map.end(); i++)
        parameters.type_map_reverse[i->second] = i->first;

    apply_parameters(systems,parameters);
    return 0;
}

void print_parameters(Parameters &parameters)
{
    for(auto &type : parameters.type_map)
    {
        printf("    Element %6s type %3d ",type.first.c_str(),type.second);
        if (parameters.lj_on || parameters.phahst_on)
            printf("eps %12.8f sig %12.8f ",parameters.lj_eps[type.second],parameters.lj_sig[type.second]);
        if (parameters.pol_on || parameters.fit_alpha_on)
            printf("alpha %12.8f ",parameters.alpha[type.second]);
        if (parameters.phahst_on)
            printf("c6 %12.8f c8 %12.7f c10 %12.6f",parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)),parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)),parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)));
        printf("\n");
    }
    printf("\n");
    return;
}

int main(int argc, char* argv[])
{
    int i;
    double error;
    vector<System> systems;
    Parameters parameters;
    printf("Force Fit\n© Adam Hogan 2021\n-------------\n");
    for (i = 0; i < argc; i++)
    {
        string arg = argv[i];
        if ( (arg == "-h") || (arg == "--help") || (argc !=2) )
        {
            printf("Usage: %s <option(s)> input\nOptions:\n-h or --help to show this dialog\n",argv[0]);
            return 0;
        }
    }

    string input_file = argv[1];
    if (load_options(systems,parameters,input_file)!=0)
    {
        printf("Error encountered loading options, exiting\n");
        return 0;
    }

    apply_parameters(systems,parameters);
    calc_forces(systems,parameters);
    printf("\nSystems loaded, beginning training ...\nInitial Error: %10.6f\nInitial parameters:\n",calc_error(systems,parameters));
    print_parameters(parameters);
    for (i=0; i<parameters.steps; i++)
    {
        bool reset = (i%20 == 0) ? true : false;
        error = train(systems, parameters, reset);
        if (i % parameters.output_freq == 0)
        {
            printf("        error: %10.6f %i/%i steps\n",error,i,parameters.steps);
            print_parameters(parameters);
        }
    }

    apply_parameters(systems,parameters);
    calc_forces(systems,parameters);
    printf("        error: %10.6f %i/%i steps\n",error,parameters.steps,parameters.steps);
    print_parameters(parameters);
    printf("\nFinal parameters\n");
    print_parameters(parameters);
    printf("\nInput parameters\n");
    for(auto &type : parameters.type_map)
    {
        printf("! %s\n",type.first.c_str());
        printf("lj_sig %d %12.8f\n",type.second,parameters.lj_sig[type.second]);
        printf("lj_eps %d %12.8f\n",type.second,parameters.lj_eps[type.second]);
        printf("c6 %d %12.8f\n",type.second,parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)));
        printf("c8 %d %12.8f\n",type.second,parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)));
        printf("c10 %d %12.8f\n",type.second,parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)));
        printf("alpha %d %12.8f\n\n",type.second,parameters.alpha[type.second]);
    }



    printf("virial parameter seds\n");
    for(auto &type : parameters.type_map)
    {
        printf("sed -i -E 's/,\"%s\",([0-9.-]*),[0-9.-]*,[0-9.-]*,[0-9.-]*,[0-9.-]*,[0-9.-]*,[0-9.-]*\\)\\)/,\"%s\",\\1,%f,%f,%f,%f,%f,%f\\)\\)/g' virial_coeff.py\n",type.first.c_str(),type.first.c_str(),parameters.alpha[type.second],parameters.lj_eps[type.second],parameters.lj_sig[type.second],parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)),parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)),parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)));
    }

    printf("\nlammps table parameters\n");
    for(auto &type : parameters.type_map)
    {
        printf("params['%s'] = {'b': %f, 'rho': %f, 'c6': %f*0.021958709/(3.166811429*0.000001), 'c8': %f*0.0061490647/(3.166811429*0.000001), 'c10': %f*0.0017219135/(3.166811429*0.000001), 'alpha': %f, 'c9': 0.0}\n",type.first.c_str(),parameters.lj_eps[type.second],parameters.lj_sig[type.second],parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)),parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)),parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)),parameters.alpha[type.second]);
    }

    printf("\nmpmc pdb seds\n");
    for(auto &type : parameters.type_map)
    {
        printf("atom type: %s --- sed -E \"s/([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)[[:space:]]*([0-9.-]*)/\\1 \\2 %8.5f %8.5f %8.5f 0.00000 0.00000 %8.5f %8.5f %8.5f/g\"\n",type.first.c_str(),parameters.alpha[type.second],parameters.lj_eps[type.second],parameters.lj_sig[type.second],parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)),parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)),parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)));
    }

    printf("\nopenff ff files\n");
    for(auto &type : parameters.type_map)
    {
        printf("%s --- sigma=\"%f * angstrom\" beta=\"%f * angstrom**-1\" c6=\"%e * kilojoule_per_mole * nanometer**6\" c8=\"%e * kilojoule_per_mole * nanometer**8\" c10=\"%e * kilojoule_per_mole * nanometer**10\">\n",type.first.c_str(),parameters.lj_sig[type.second],parameters.lj_eps[type.second],parameters.c6[type.second] / (0.021958709 / (3.166811429 * 0.000001)) * 1.0/0.000003166811563*0.008314462618/pow(10*1.88973,6),parameters.c8[type.second] / (0.0061490647 / (3.166811429 * 0.000001)) * 1.0/0.000003166811563*0.008314462618/pow(10*1.88973,8),parameters.c10[type.second] / (0.0017219135 / (3.166811429 * 0.000001)) * 1.0/0.000003166811563*0.008314462618/pow(10*1.88973,10));
    }
    for(auto &type : parameters.type_map)
    {
        printf("%s --- alpha=\"%f * angstrom**3\">\n",type.first.c_str(), parameters.alpha[type.second]);
    }

    printf("\n\nFinal fit\n\n%6s --- %25s %25s %25s\n","Config","ab initio E","Model E","Error");
    if (parameters.energy_fitting_on)
    {
        int i = 0;
        for (auto &system : systems)
        {
            double fit_energy = system.fit_energy <= 0. ? system.fit_energy : parameters.max_energy * atan(system.fit_energy/parameters.max_energy);
            double calc_energy = system.total_energy <= 0. ? system.total_energy : parameters.max_energy * atan(system.total_energy/parameters.max_energy);
            double error = pow(calc_energy-fit_energy,2);
            printf("%6d --- %25.6f %25.6f %25.6f\n",++i,system.fit_energy,system.total_energy,error);
        }
    }
    
    return 0;
}

