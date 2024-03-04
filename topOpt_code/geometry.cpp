#include "geometry.h"

//--------------------------------------------
// PHYSICS CLASS
//--------------------------------------------
void PHYSICS::initialize()
{
    eval_solution_times();
    build_bounds_elems_v();
    build_bounds_elems_surfaces_v();
    
    // for (int ibound = 0; ibound < nBounds; ibound++)
    // {
    //     bounds_elems_v[ibound].print();
    //     bounds_elems_surface_v[ibound].printRowMatlab(std::to_string(ibound));
    // }
}

void PHYSICS::build_bounds_elems_v()
{
    bounds_elems_v.resize(nBounds);
    std::string folder_path = name;
    folder_path = "PREPRO/PROBLEM_DATA/" + folder_path;
    std::string file_path = folder_path + "/BoundElems_v.txt";
    std::ifstream  bound_file; 
    bound_file.open(file_path); 
    if (!bound_file.is_open()) 
    {
        throw_line("ERROR, can't open input data file");
    }
    std::string line;
    std::istringstream iss;
    
    for (int ibound = 0; ibound < nBounds; ibound++)
    {
        STREAM::getLines(bound_file, line, 1);
        VECTOR_INT check(2);
        STREAM::getRowVector(bound_file, line, iss, check);
        if (check[0] != ibound)
        {
            throw_line("ERROR: uncompatible ibound and reading bound\n");
        }
        int n_el_in_bound = check[1];
        bounds_elems_v[ibound].initialize(n_el_in_bound, dim);
        STREAM::getMatrix(bound_file, line, iss, bounds_elems_v[ibound]);
    }
    bound_file.close();
}

void PHYSICS::build_bounds_elems_surfaces_v()
{
    bounds_elems_surface_v.resize(nBounds);
    for (int ibound = 0; ibound < nBounds; ibound++)
    {
        MATRIX_INT bound_elems = bounds_elems_v[ibound];
        bounds_elems_surface_v[ibound].initialize(bound_elems.nRow);
        for (int iel = 0; iel < bound_elems.nRow; iel++)
        {
            MATRIX coord(dim, dim);
            for (int inode = 0; inode < dim; inode++)
            {
                int node = bound_elems[iel][inode];
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    coord[inode][jcomp] = coord_v[node][jcomp];
                }
            }
            bounds_elems_surface_v[ibound][iel] = get_surface(coord, dim);
        }
    }
}

void PHYSICS::print()
{
    printf(" dim: %d\n rho: %" format "\n: mu: %" format" \n ni: %" format "\n nNodes_v:%d\n", dim, rho, mu, ni, nNodes_v);
}

void PHYSICS::eval_solution_times()
{
    switch (isStationary)
    {
        case 1: // stationary solution
        {
            t_end = 1e16;
            deltaT = t_end;
            deltaT_min = 1e-5;
            solution_times.append(t_end);
            solution_deltaT.append(deltaT);
            break;
        }
        default: // time dependent solution
        {
            solution_times.resetZeros();
            prec curr_time = 0.0;
            if (abs(t_end-curr_time) <= 1e-10) throw_line("ERROR: too small ending time (<= 1e-10)\n");
            while (abs(t_end-curr_time) > 1e-10)
            {
                solution_times.append(curr_time);
                curr_time += deltaT;
            }
            solution_times.append(t_end);
            solution_deltaT.initialize(solution_times.length-1);
            solution_deltaT.reset(deltaT);
            solution_deltaT[solution_deltaT.length-1] = t_end - solution_times[solution_times.length-2]; //adjust the last time step
            break;
        }
    }
}
//--------------------------------------------
// END PHYSICS CLASS
//--------------------------------------------

//--------------------------------------------
// BASE CONSTRAINT CLASS
//--------------------------------------------
void CONSTRAINT::initialize(int constraint_type)
{
    if (constraint_type < 0)
    {
        throw_line ("ERROR: invalid constraint type\n");
    }
    else
    {
        type = constraint_type;
    }
}

void CONSTRAINT::initialize(int constraint_type, VECTOR constraint_parameters)
{
    int max_type = 0;
    if ((constraint_type < 0) || (constraint_type > max_type))
    {
        throw_line ("ERROR: invalid constraint type\n");
    }
    else
    {
        type = constraint_type;
    }
    parameters = constraint_parameters;
}

void CONSTRAINT::initialize_volume_constraint(int constr_type, prec volume_fraction)
{
    initialize(constr_type);
    Vr = volume_fraction;
}

void CONSTRAINT::initialize_subdomain_volume_constraint(int constr_type, int subdomain, prec volume_fraction)
{
    initialize_volume_constraint(constr_type, volume_fraction);
    domain_id = subdomain;
}

void CONSTRAINT::initialize_bound_constraint(int constr_type, int bound, prec bound_fraction)
{
    initialize(constr_type);
    bound_id = bound;
    Sr = bound_fraction;
}
//--------------------------------------------
// END BASE CONSTRAINT CLASS
//--------------------------------------------

//--------------------------------------------
// CONSTRAINTS HANDLER CLASS
//--------------------------------------------
void CONSTRAINTS::initialize(int constraints_n_constr, VECTOR_INT constraints_types, std::vector<VECTOR> constraints_parameters)
{
    if (constraints_n_constr < 1)
    {
        throw_line("ERROR: invalid number of constraints. Less than 1\n");
    }
    n_constr = constraints_n_constr;
    save_constraints_type(constraints_types);
    build_constraints_list(constraints_parameters);
}

void CONSTRAINTS::save_constraints_type(VECTOR_INT constraints_types)
{
    types.setZeros(n_constr);
    if (constraints_types.length != n_constr)
    {
        std::cout << "\nn_constr: " << n_constr << "\t types_length: " << constraints_types.length << "\n";
        throw_line("ERROR: incoherent number of contraints and constraint types\n");
    }
    for (int icons = 0; icons < n_constr; icons++)
    {
        int constraint_type = constraints_types[icons];
        if (( constraint_type < 0) || (constraint_type > max_constraint_type))
        {
            std::cout << "\nconstraint: " << icons+1 << " / " << n_constr << "\n";
            std::cout << "type: " << constraint_type << "\n";
            throw_line("ERROR: not implemented constraint type\n");
        }
        else
        {
            types[icons] = constraint_type;
        }
    }
}

void CONSTRAINTS::build_constraints_list(std::vector<VECTOR> constraints_parameters)
{
    if (int(constraints_parameters.size()) != n_constr)
    {
        std::cout << "\nn_constr: " << n_constr << "\t parameters_length: " << constraints_parameters.size() << "\n";
        throw_line("ERROR: incoherent number of contraints and constraint types\n");
    }
    list.resize(n_constr);
    for (int icons = 0; icons < n_constr; icons++)
    {
        int type = types[icons];
        switch (type)
        {
            case 0:
            {
                prec volume_fraction = constraints_parameters[icons][0];
                CONSTRAINT vol_constraint;
                vol_constraint.initialize_volume_constraint(type, volume_fraction);
                list[icons] = vol_constraint;
                break;
            }  
            case 1:
            {
                int subdomain = int(constraints_parameters[icons][0]);
                prec volume_fraction = constraints_parameters[icons][1];
                CONSTRAINT subdom_vol_constraint;
                subdom_vol_constraint.initialize_subdomain_volume_constraint(type, subdomain, volume_fraction);
                list[icons] = subdom_vol_constraint;
                break;
            } 
            case 2:
            {
                int bound = int(constraints_parameters[icons][0]);
                prec bound_fraction = constraints_parameters[icons][1];
                CONSTRAINT bound_constraint;
                bound_constraint.initialize_bound_constraint(type, bound, bound_fraction);
                list[icons] = bound_constraint;
                break;
            }
            default:
            {
                throw_line("ERROR: unexpected invalid constraint type, this error should have already popUp in types initialzation\n");
                break;
            }
        }
    }
}
//--------------------------------------------
// END CONSTRAINTS HANDLER CLASS
//--------------------------------------------
