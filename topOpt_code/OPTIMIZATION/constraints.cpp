#include "constraints.h"

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

void CONSTRAINT::initialize_discretizing_constraint(int constr_type, prec discretization_tollerance)
{
    initialize(constr_type);
    discretization_toll = discretization_tollerance;
}

void CONSTRAINT::initialize_WSS_constraint(int constr_type, prec crit_WSS, int sign_value)
{
    // sign = -1: lower than crit_WSS
    // sign = 1: greater than crit_WSS
    initialize(constr_type);
    critical_WSS = crit_WSS;
    if ((sign_value != -1) && (sign_value != 1))
    {
        std::cout << "sign: " << sign_value << "\n";
        throw_line("ERROR: not valid sign value initializing WSS constraint\n");
    }
    sign = sign_value;
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
            case 3:
            {
                prec discretization_toll = constraints_parameters[icons][0];
                CONSTRAINT discretizing_constraint;
                discretizing_constraint.initialize_discretizing_constraint(type, discretization_toll);
                list[icons] = discretizing_constraint;
                break;
            }
            case 4:
            {
                prec crit_WSS = constraints_parameters[icons][0];
                CONSTRAINT WSS_constraint;
                int sign_value = int(constraints_parameters[icons][1]);
                WSS_constraint.initialize_WSS_constraint(type, crit_WSS, sign_value);
                list[icons] = WSS_constraint;
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

//--------------------------------------------
// START WSS CONSTRAINT HANDLER CLASS
//--------------------------------------------
void WSS_CONSTRAINT::initialize(PHYSICS &physics_ptr, std::vector<MATRIX> velocity_gradient_in, std::vector<VECTOR> gamma_gradient_in, std::vector<std::vector<VECTOR>> gamma_tangents_in)
{
    physics = &physics_ptr;
    velocity_gradient = velocity_gradient_in;
    gamma_gradient = gamma_gradient_in;
    gamma_tangents = gamma_tangents_in;
    nodal_WSS_tangent_values.resize((*physics).nNodes_v);
    n_tan =  int(gamma_tangents[0].size());
    if (n_tan != ((*physics).dim-1))
    {
        throw_line("ERROR: tagnent variety dimension != (dim-1)");
    }
    for (int inod = 0; inod < (*physics).nNodes_v; inod++)
    {
        nodal_WSS_tangent_values[inod].setZeros(n_tan);
    }
    nodal_WSS_values.setZeros((*physics).nNodes_v);
}

void WSS_CONSTRAINT::eval_nodal_values()
{
    int nNodes_v = (*physics).nNodes_v;
    prec mu = (*physics).mu;
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        VECTOR temp_vec = velocity_gradient[inod] * gamma_gradient[inod];
        for (int itan = 0 ; itan < n_tan; itan++)
        {
            prec temp_val = mu * temp_vec.dot(gamma_tangents[inod][itan]); 
            nodal_WSS_tangent_values[inod][itan] = temp_val; // = mu * [(grad(U)_i * grad(gamma)_i)^t * grad_gamma_tangent_i]
            nodal_WSS_values[inod] += temp_val * temp_val;
        }
        nodal_WSS_values[inod] = sqrt(nodal_WSS_values[inod]); // WSS = sum_over_tangents("tangent_WSS" ^ 2)
    }
}
//--------------------------------------------
// END WSS CONSTRAINT HANDLER CLASS
//--------------------------------------------