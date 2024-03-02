#include "geometry.h"

//--------------------------------------------
// PHYSICS CLASS
//--------------------------------------------
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
    int max_type = 0;
    if ((constraint_type < 0) || (constraint_type > max_type))
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

void CONSTRAINT::initialize_edge_constraint(int constr_type, int edge, prec edge_fraction)
{
    initialize(constr_type);
    edge_id = edge;
    Sr = edge_fraction;
}
//--------------------------------------------
// END BASE CONSTRAINT CLASS
//--------------------------------------------

// //--------------------------------------------
// // VOLUME CONSTRAINT CLASS
// //--------------------------------------------
// void VOLUME_CONSTRAINT::initialize(int constraint_type, VECTOR constraint_parameters)
// {
//     // initialize constraint
//     CONSTRAINT::initialize(constraint_type);

//     // define maximum volume fraction
//     prec constraint_Vr = constraint_parameters[0];
//     if ((0.0 >= constraint_Vr) || (constraint_Vr > 1.0))
//     {
//         std::cout << "\nVr: " << constraint_Vr << "\n";
//         throw_line ("ERROR: invalid maximum volume fraction\n");
//     }
//     else
//     {
//         Vr = constraint_Vr; // the Vr  must be inserted in the first index
//     }
// }
// //--------------------------------------------
// // END VOLUME CONSTRAINT CLASS
// //--------------------------------------------

// //--------------------------------------------
// // EDGE SIZE CONSTRAINT CLASS
// //--------------------------------------------
// void EDGE_SIZE_CONSTRAINT::initialize(int constraint_type, VECTOR constraint_parameters)
// {
//     // initialize constraint
//     CONSTRAINT::initialize(constraint_type);

//     // define edge_id
//     int constraint_edge_id = int(constraint_parameters[0]);
//     edge_id = constraint_edge_id;

//     // define maximum edge size fraction
//     prec constraint_Sr = constraint_parameters[1];
//     if ((0.0 > constraint_Sr) || (constraint_Sr > 1.0))
//     {
//         throw_line ("ERROR: invalid maximum volume fraction\n");
//     }
//     else
//     {
//         Sr = constraint_Sr; // the Vr  must be inserted in the first index
//     }
// }
// //--------------------------------------------
// // END EDGE SIZE CONSTRAINT CLASS
// //--------------------------------------------

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
            throw_line("ERROR: invalid constraint type\n");
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
                    int edge = int(constraints_parameters[icons][0]);
                    prec edge_fraction = constraints_parameters[icons][1];
                    CONSTRAINT edge_constraint;
                    edge_constraint.initialize_edge_constraint(type, edge, edge_fraction);
                    list[icons] = edge_constraint;
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
