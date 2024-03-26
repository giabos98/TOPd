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

    build_centroids_v();
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

//build elements centorids
void PHYSICS::build_centroids_v()
{
    centroids_v.resize(nElem_v);
    for (int iel = 0; iel < nElem_v; iel++)
    {
        centroids_v[iel].initialize(dim);
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            VECTOR temp_coord = coord_v.get_row(iglob);
            centroids_v[iel] += temp_coord;
        }
        prec n_nodes_x_el = (dim+1) * 1.0;
        centroids_v[iel] /= n_nodes_x_el;
    }
}

//eval a quantity on elements
void PHYSICS::eval_on_elements_v(VECTOR &value, VECTOR &element_value)
{
    // with 1D schape functions, the average value on the element corresponds to the aritmetic mean of the nodal values
    // so it is fine to sum the nodal value for each element, and finally divide th ewhole vector by th number of nodes per element
    // value is supposed known in the coord_v mesh
    element_value.completeReset();
    element_value.setZeros(nElem_v);
    for (int iel = 0; iel < nElem_v; iel++)
    {
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            element_value[iel] += value[iglob];
        }
    }
    prec n_nodes_x_el = (dim+1) * 1.0; // make dim+1 a prec variable to divide it to element_value
    element_value /= n_nodes_x_el;
}

//eval a quantity on elements centroids
void PHYSICS::eval_on_centroids_v(VECTOR &value, VECTOR &centroids_value)
{
    // value is supposed known in the coord_v mesh
    centroids_value.completeReset();
    centroids_value.setZeros(nElem_v);
    for (int iel = 0; iel < nElem_v; iel++)
    {
        prec temp_weights_sum = 0.0;
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            VECTOR temp_coords = coord_v.get_row(iglob);
            prec temp_dist = (centroids_v[iel] - temp_coords).norm();
            temp_weights_sum += temp_dist;
            centroids_value[iel] += value[iglob] * temp_dist;
        }
        centroids_value[iel] /= temp_weights_sum;
    }
}

// static prec get_surface(MATRIX &matCoord, int dim);
prec PHYSICS::get_surface(MATRIX &matCoord, int dim)
{
    if (dim != 2 && dim != 3) throw_line("ERROR: Getting area in a space different from 2D or 3D");
    if (matCoord.nRow != dim || matCoord.nCol != dim) throw_line("ERROR: Getting area of coords incompatible with dim\n");
    prec area;
    switch (dim)
    {
        case 2:
        {
            std::shared_ptr<prec[]>  vec(new prec[dim]);
            vec[0] = matCoord[0][0] - matCoord[1][0];
            vec[1] = matCoord[0][1] - matCoord[1][1];
            area = VECTOR::norm(vec,2);
            break;
        }
        case 3:
        {
            VECTOR vec31(dim);
            vec31[0] = matCoord[2][0] - matCoord[0][0];
            vec31[1] = matCoord[2][1] - matCoord[0][1];
            vec31[2] = matCoord[2][2] - matCoord[0][2];
            VECTOR vec21(dim);
            vec21[0] = matCoord[1][0] - matCoord[0][0];
            vec21[1] = matCoord[1][1] - matCoord[0][1];
            vec21[2] = matCoord[1][2] - matCoord[0][2];
            // VECTOR tempVec; tempVec = vec31.cross(vec21);
            // vec21 = tempVec;
            VECTOR newVec;
            VECTOR::cross(vec31, vec21, newVec);
            
            area = VECTOR::norm(newVec)/2;
            break;
        }
    }
    return area;
}

void PHYSICS::eval_gradient(VECTOR &value, std::vector<VECTOR> &gradient, int eval_method)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[nodal_value]
    // gradient structure: gradient[derivative_component][nodel_value]
    //eval_method = 0: Lampeato. Other cases not implemented
    if (eval_method == 0) // Lampeato
    {
        gradient.resize(dim);
        VECTOR weights; 
        weights.setZeros(nNodes_v);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            gradient[icomp].completeReset();
            gradient[icomp].setZeros(value.length);
        }
        for (int iel = 0; iel < nElem_v; iel++)
        {
            prec temp_weight = Volume_v[iel];
            prec element_value = 0.0;
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[iel][iloc];
                weights[iglob] += temp_weight;
                element_value += value[iglob];
            }
            element_value /= (dim+1)*1.0;
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[iel][iloc];
                prec nodal_value = value[iglob];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    if ((iglob == 2392))
                    {
                        std::cout << "\nicomp: " << icomp << "\n";
                        std::cout << "iglob: " << iglob << "\n";
                        std::cout << "pre-grad: " << gradient[icomp][iglob] << "\n";
                        std::cout << "nodal: " << nodal_value << "\n";
                        std::cout << "elem: " << element_value << "\n";
                        std::cout << "coef: " << Coef_v[icomp][iel][iloc] << "\n";
                        std::cout << "weight: " << temp_weight << "\n";
                        std::cout << "INCREMENT: " << element_value * Coef_v[icomp][iel][iloc] * temp_weight << "\n";
                    }
                    gradient[icomp][iglob] += element_value * temp_weight * Coef_v[icomp][iel][iloc];
                    if ((iglob == 2392))
                    {
                        std::cout << "post-grad: " << gradient[icomp][iglob] << "\n";
                        // pause();

                    }
                }
            }
        }
        std::cout << "mid_grad: " << gradient[0][2392] << "\n";
        // pause();
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int iglob = 0; iglob < nNodes_v; iglob++)
            {
                if ((iglob == 2392))
                {
                    std::cout << "\npre-comp " << icomp << ": " << gradient[icomp][iglob] << "\n";
                }
                gradient[icomp][iglob] /= weights[iglob];
                if ((iglob == 2392))
                {
                    std::cout << "\npost-comp " << icomp << ": " << gradient[icomp][iglob] << "\n";
                    // pause();
                }
            }
        }
    }
    else
    {
        throw_line("ERROR: not handled gradient reconstruction method\n");
    }
}

void PHYSICS::eval_gradient(MATRIX &value, std::vector<std::vector<VECTOR>> &gradient, int eval_method)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[component][nodal_value]
    // gradient structure: gradient[value_component][derivative_component][nodel_value]
    //initialize solution
    gradient.resize(dim);
    int n_comps = value.nRow;
    for (int icomp = 0; icomp < n_comps; icomp++)
    {
        VECTOR comp_value = value.get_row(icomp);
        eval_gradient(comp_value, gradient[icomp], eval_method);
    }
}

void PHYSICS::eval_gradient_norm(std::vector<VECTOR> &gradient, VECTOR &norm)
{
    norm.completeReset();
    norm.setZeros(nNodes_v);
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        for (int icomp = 0; icomp < dim; icomp++)
        {
            norm[inod] += gradient[icomp][inod]*gradient[icomp][inod];
        }
    }
    norm.squared_root();
}

void PHYSICS::eval_gradient_norm(std::vector<std::vector<VECTOR>> &gradient, VECTOR &norm)
{
    norm.completeReset();
    norm.setZeros(nNodes_v);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        std::vector<VECTOR> comp_gradient = gradient[icomp];
        VECTOR comp_norm(nNodes_v);
        eval_gradient_norm(comp_gradient, comp_norm);
        comp_norm.power(2);
        norm += comp_norm;
    }
    norm.squared_root();
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

void CONSTRAINT::initialize_discretizing_constraint(int constr_type, prec discretization_tollerance)
{
    initialize(constr_type);
    discretization_toll = discretization_tollerance;
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
                prec discretization_toll = int(constraints_parameters[icons][0]);
                CONSTRAINT discretizing_constraint;
                discretizing_constraint.initialize_discretizing_constraint(type, discretization_toll);
                list[icons] = discretizing_constraint;
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
