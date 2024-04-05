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
    parse_bounds();
    build_centroids_v();
}

void PHYSICS::parse_bounds()
{
    std::string folderPath = name;
    folderPath = "PREPRO/PROBLEM_DATA/" + folderPath;
    bound_elems_v_path = folderPath + "/BoundElems_v.txt";
    bound_nodes_v_path = folderPath + "/BoundNodes_v.txt";

    // BEGIN NODES_V STREAMING
    std::ifstream bound_nodes_v_stream;
    bound_nodes_v_stream.open(bound_nodes_v_path);
    if (!bound_nodes_v_stream.is_open()) throw_line("ERROR, can't open input data file");
    std::string line;
    std::istringstream iss;

    VECTOR_INT general_info(2);
    STREAM::getRowVector(bound_nodes_v_stream, line, iss, general_info);
    int n_bounds = general_info[1];
    bound_nodes_v.resize(n_bounds);

    for (int ibound = 0; ibound < n_bounds; ibound++)
    {
        VECTOR_INT temp_info(2);
        STREAM::getRowVector(bound_nodes_v_stream, line, iss, temp_info);
        int temp_n_nodes = temp_info[1];
        bound_nodes_v[ibound].initialize(temp_n_nodes);
        VECTOR_INT temp_bound_nodes(temp_n_nodes);
        STREAM::getColVector(bound_nodes_v_stream, line, iss, temp_bound_nodes, temp_bound_nodes.length);
        bound_nodes_v[ibound] = temp_bound_nodes;
        getline(bound_nodes_v_stream, line);
    }  
    bound_nodes_v_stream.close();
    // CLOSE NODES_V STREAMING

    // BEGIN ELEMS_V STREAMING
    std::ifstream bound_elems_v_stream;
    bound_elems_v_stream.open(bound_elems_v_path);
    if (!bound_elems_v_stream.is_open()) throw_line("ERROR, can't open input data file");
    // std::string line;
    // std::istringstream iss;

    getline(bound_elems_v_stream, line); // skip general info
    bound_elems_v.resize(n_bounds);

    for (int ibound = 0; ibound < n_bounds; ibound++)
    {
        VECTOR_INT temp_info(2);
        STREAM::getRowVector(bound_elems_v_stream, line, iss, temp_info);
        int temp_n_elem = temp_info[1];
        bound_elems_v[ibound].initialize(temp_n_elem, dim);
        MATRIX_INT temp_bound_elems(temp_n_elem, dim);
        STREAM::getMatrix(bound_elems_v_stream, line, iss, temp_bound_elems);
        bound_elems_v[ibound] = temp_bound_elems;
        getline(bound_elems_v_stream, line);
    } 
    bound_elems_v_stream.close();
    // CLOSE ELEMS_V STREAMING   
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
            t_end = 1e10;
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

//eval a quantity on the elements of the velocity mesh
void PHYSICS::eval_on_elements_v(VECTOR &value, VECTOR &element_value)
{
    // with 1D shape functions, the average value on the element corresponds to the aritmetic mean of the nodal values
    // so it is fine to sum the nodal value for each element, and finally divide the whole vector by th number of nodes per element
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

//eval a quantity on a given subset 
void PHYSICS::eval_on_elements_v(VECTOR &value, VECTOR_INT &nodesFromNodes_v, VECTOR_INT &elems, VECTOR &element_value)
{
    // with 1D shape functions, the average value on the element corresponds to the aritmetic mean of the nodal values
    // so it is fine to sum the nodal value for each element, and finally divide the whole vector by th number of nodes per element
    // value is supposed known in the coord_v mesh nodes belonging to the given elements, therefore its length will be the number of different nodes in the given subset of elements
    // nodesFromNodes_v is used to recall the 'id' of a subset node from its global one
    element_value.completeReset();
    element_value.setZeros(elems.length);
    for (int iel = 0; iel < elems.length; iel++)
    {
        int glob_el = elems[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int glob_id = elem_v[glob_el][iloc];
            int sub_glob_id = nodesFromNodes_v[glob_id];
            element_value[iel] += value[sub_glob_id];
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

// smooth a value between the elements of the velocity mesh
void PHYSICS::smooth_between_elements(VECTOR &value, VECTOR &smoothed_value)
{
    // value is supposed know in all the velocity mesh nodes
    smoothed_value.completeReset();
    smoothed_value.setZeros(nNodes_v);

    VECTOR value_on_elems;
    value_on_elems.setZeros(nElem_v);
    eval_on_elements_v(value, value_on_elems);

    VECTOR nodal_weights;
    nodal_weights.setZeros(nNodes_v);
    for (int iel = 0; iel < nElem_v; iel++)
    {
        prec temp_volume = Volume_v[iel];
        prec temp_elem_value = value_on_elems[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            smoothed_value[iglob] += temp_elem_value * temp_volume;
            nodal_weights[iglob]  += temp_volume;
        }
    }
    smoothed_value *= nodal_weights;
}

// smooth a value between a given subset of elements of the velocity mesh
void PHYSICS::smooth_between_elements(VECTOR &value, VECTOR_INT &nodesFromNodes_v, VECTOR_INT &elems, VECTOR &smoothed_value)
{
    // value is supposed known in the coord_v mesh nodes belonging to the given elements, therefore its length will be the number of different nodes in the given subset of elements
    // nodesFromNodes_v is used to recall the 'id' of a subset node from its global one
    smoothed_value.completeReset();
    smoothed_value.setZeros(value.length);

    VECTOR value_on_elems;
    value_on_elems.setZeros(elems.length);
    eval_on_elements_v(value, nodesFromNodes_v, elems, value_on_elems);
    // value_on_elems.printRowMatlab("elems");

    VECTOR nodal_weights;
    nodal_weights.setZeros(value.length);
    for (int iel = 0; iel < elems.length; iel++)
    {
        int glob_el = elems[iel];
        prec temp_volume = Volume_v[glob_el];
        prec temp_elem_value = value_on_elems[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int glob_id = elem_v[glob_el][iloc];
            int sub_glob_id = nodesFromNodes_v[glob_id];
            smoothed_value[sub_glob_id] += temp_elem_value * temp_volume;
            nodal_weights[sub_glob_id]  += temp_volume;
        }
    }
    smoothed_value /= nodal_weights;
    // smoothed_value.printRowMatlab("smooth");
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
    // gradient structure: gradient[derivative_component][nodal_value]
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
            VECTOR element_gradient_value; element_gradient_value.setZeros(dim);
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[iel][iloc];
                weights[iglob] += temp_weight;
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    element_gradient_value[icomp] += value[iglob]*Coef_v[icomp][iel][iloc];
                }
            }
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[iel][iloc];
                // prec nodal_value = value[iglob];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    // if ((iglob == 2392))
                    // {
                    //     std::cout << "\nicomp: " << icomp << "\n";
                    //     std::cout << "iglob: " << iglob << "\n";
                    //     std::cout << "pre-grad: " << gradient[icomp][iglob] << "\n";
                    //     std::cout << "nodal: " << nodal_value << "\n";
                    //     std::cout << "elem: " << element_value << "\n";
                    //     std::cout << "coef: " << Coef_v[icomp][iel][iloc] << "\n";
                    //     std::cout << "weight: " << temp_weight << "\n";
                    //     std::cout << "INCREMENT: " << element_value * Coef_v[icomp][iel][iloc] * temp_weight << "\n";
                    // }
                    gradient[icomp][iglob] += element_gradient_value[icomp] * temp_weight;
                    // if ((iglob == 2392))
                    // {
                    //     std::cout << "post-grad: " << gradient[icomp][iglob] << "\n";
                    //     // pause();

                    // }
                }
            }
        }
        // std::cout << "mid_grad: " << gradient[0][2392] << "\n";
        // pause();
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int iglob = 0; iglob < nNodes_v; iglob++)
            {
                // if ((iglob == 2392))
                // {
                //     std::cout << "\npre-comp " << icomp << ": " << gradient[icomp][iglob] << "\n";
                // }
                gradient[icomp][iglob] /= weights[iglob];
                // if ((iglob == 2392))
                // {
                //     std::cout << "\npost-comp " << icomp << ": " << gradient[icomp][iglob] << "\n";
                //     // pause();
                // }
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
    // gradient structure: gradient[value_component][derivative_component][nodal_value]
    //initialize solution
    gradient.resize(dim);
    int n_comps = value.nRow;
    for (int icomp = 0; icomp < n_comps; icomp++)
    {
        VECTOR comp_value = value.get_row(icomp);
        eval_gradient(comp_value, gradient[icomp], eval_method);
    }
}

void PHYSICS::eval_directional_gradient(VECTOR &value, VECTOR_INT &nodes, std::vector<VECTOR> &directions, VECTOR &dir_gradient)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[nodal_value]
    // nodes structure: nodes[specific_node] = global_node
    // directions structure: directions[specific_node][direction]
    // dir_gradient structure: gradient[specific_nodal_value]
    // N.B.: specific node means it belongs to [0,nodes-length-1]
    dir_gradient.completeReset();
    dir_gradient.setZeros(nodes.length);
    std::vector<VECTOR> gradient;
    eval_gradient(value, gradient);
    for (int inod = 0; inod < nodes.length; inod++)
    {
        int iglob = nodes[inod];
        VECTOR temp_dir = directions[inod];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            dir_gradient[inod] += gradient[icomp][iglob] * temp_dir[icomp];
        }
    }
}

void PHYSICS::eval_directional_gradient(MATRIX &value, VECTOR_INT &nodes, std::vector<VECTOR> &directions, std::vector<VECTOR> &dir_gradient)
{  
    // suppose the value interpolated in the coord_v
    // value structure: value[component][nodal_value]
    // nodes structure: nodes[specific_node] = global_node
    // directions structure: directions[specific_node][direction]
    // dir_gradient structure: gradient[value_component][specific_nodal_value]
    // N.B.: specific node means it belongs to [0,nodes-length-1]

    dir_gradient.resize(dim);
    int n_comps = value.nRow;
    for (int icomp = 0; icomp < n_comps; icomp++)
    {
        VECTOR comp_value = value.get_row(icomp);
        eval_directional_gradient(comp_value, nodes, directions, dir_gradient[icomp]);
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

void PHYSICS::eval_WSS(MATRIX &value, VECTOR_INT &nodes, std::vector<VECTOR> &normals, VECTOR &WSS)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[component][nodal_value]
    // normals structure: normals[node][component]
    // WSS structure: WSS[specific_nodal_value]
    
    if (nodes.length != int(normals.size()))
    {
        throw_line("ERROR: the nodes and the normals for the WSS evaluation are defined in a different number of points\n");
    }
    int n_nodes = nodes.length;
    WSS.completeReset(); 
    WSS.setZeros(nodes.length);
    std::vector<VECTOR> normal_gradient;
    eval_directional_gradient(value, nodes, normals, normal_gradient);
    int n_tan = dim-1;
    std::vector<std::vector<VECTOR>> tangent_vectors(n_nodes);

    // build tangent vectors
    for (int inod = 0; inod < n_nodes; inod++)
    {
        tangent_vectors[inod].resize(n_tan);
        VECTOR temp_normal = normals[inod];
        switch (n_tan)
        {
            case 1: // 2D
            {
                tangent_vectors[inod][0].setZeros(dim);
                tangent_vectors[inod][0][0] = -temp_normal[1];
                tangent_vectors[inod][0][1] = temp_normal[0];
                tangent_vectors[inod][0] /= tangent_vectors[inod][0].norm();
                break;
            }
            case 2: // 3D
            {
                // first tangent vector
                tangent_vectors[inod][0].setZeros(dim);
                tangent_vectors[inod][0][0] = -temp_normal[1];
                tangent_vectors[inod][0][1] = temp_normal[0];
                prec first_tan_norm = tangent_vectors[inod][0].norm();
                tangent_vectors[inod][0] /= first_tan_norm;

                // second tangent vector, orthogonal to the first one
                tangent_vectors[inod][1].setZeros(dim);
                tangent_vectors[inod][1][0] = -temp_normal[2] + (temp_normal[1]*temp_normal[1]*temp_normal[2])/(first_tan_norm*first_tan_norm);
                tangent_vectors[inod][1][1] = -(temp_normal[0]*temp_normal[1]*temp_normal[2])/(first_tan_norm*first_tan_norm);
                tangent_vectors[inod][1][2] = temp_normal[0];
                tangent_vectors[inod][1] /= tangent_vectors[inod][1].norm();
                break;
            }   
            default:
            {
                throw_line("ERROR with dimensions\n");
                break;
            }
        }
    }

    // project normal gradient into tangent direction(2D)/plane(3D)
    for (int inod = 0; inod < n_nodes; inod++)
    {
        prec temp_shear_rate_norm = 0;
        VECTOR nodal_grad(dim);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            nodal_grad[icomp] = normal_gradient[icomp][inod];
        } 
        for (int itan = 0; itan < n_tan; itan++)
        {
            prec temp_grad_proj = nodal_grad.dot(tangent_vectors[inod][itan]);
            temp_shear_rate_norm += temp_grad_proj*temp_grad_proj;
        }
        temp_shear_rate_norm = sqrt(temp_shear_rate_norm);
        WSS[inod] = mu * temp_shear_rate_norm;
    }
}

void PHYSICS::eval_WSS(MATRIX &value, int bound_id, VECTOR normal, VECTOR &WSS)
{
    // the function evaluates the WSS in the nodes belonging to a given bound_id, considering as the fixed normal to each of the nodes the provided one
    // WSS structure : WSS[specific_nodal_value_in_bound_id]
    VECTOR_INT nodes = bound_nodes_v[bound_id]; 
    // !! MISSING NODES RECONSTRUCTION !!
    std::vector<VECTOR> normals(nodes.length);
    for (int inod = 0; inod < nodes.length; inod++)
    {
        normals[inod] = normal;
    }
    eval_WSS(value, nodes, normals, WSS);
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
