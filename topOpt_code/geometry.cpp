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

void PHYSICS::set_alpha_min(VECTOR &alpha_min_vector) 
{
    prec a_min = alpha_min_vector[0];
    if (dim == 2) // for the 3D case it must be set equal to 0, tha tis the default value
    {
        half_domain_thickness = alpha_min_vector[1];
    }
    prec h = half_domain_thickness;

    if (h < 0.0)
    {
        throw_line("ERROR: invalid half out of plane domain thickness\n");
    }
    else
    {
        if (h == 0)
        {
            alpha_min = a_min;
        }
        else
        {
            prec a_theory = 5*mu / (2*h*h); // folloving Borrval & Perersson (2003), Gersborg-Hansen et al. (2005)
            if (a_min >= a_theory)
            {
                alpha_min = a_min;
            }
            else
            {
                alpha_min = a_theory; 
            }
        }  
    }
}

void PHYSICS::update_real_alpha_max()
{
    // scale alpha max for the given number of iteration. USEFUL TO ALLOW GENERAL GEOMETRY FORMATION
    if (curr_opt_it < alpha_it)
    {
        real_alpha_max = ((1.0 * curr_opt_it) / (1.0 * alpha_it) * (alpha_max - alpha_min)) + alpha_min;
    }
    else
    {
        real_alpha_max = alpha_max;
    }
}

void PHYSICS::update_alpha(VECTOR &gamma, VECTOR &alpha)
{
    update_real_alpha_max();
    prec factor = q*(real_alpha_max - alpha_min);
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        alpha[inod] = alpha_min + factor * (1 - gamma[inod]) / ( q + gamma[inod]);
    }
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

void PHYSICS::decompose_NS_solution(VECTOR &sol, MATRIX &U_sol, VECTOR &P_sol, int transpose_U)
{
    U_sol.complete_reset(); P_sol.complete_reset();
    U_sol.setZeros(dim, nNodes_v);
    P_sol.setZeros(nNodes);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        int start_v_comp_id = icomp * nNodes_v;
        for (int inode = 0; inode < nNodes_v; inode++)
        {
            U_sol[icomp][inode] = sol[start_v_comp_id + inode];
        }
    }
    int start_p_sol_id = dim*nNodes_v;
    for (int inode = 0; inode < nNodes; inode++)
    {
        P_sol[inode] = sol[start_p_sol_id + inode];
    }
}

void PHYSICS::decompose_NS_solution_over_time(MATRIX &sol, std::vector<MATRIX> &U_sol, std::vector<VECTOR> &P_sol)
{
    int n_times = sol.nRow;
    U_sol.clear(); P_sol.clear();
    U_sol.resize(n_times); P_sol.resize(n_times);
    for (int itime = 0; itime < n_times; itime++)
    {
        VECTOR temp_sol = sol.get_row(itime);
        decompose_NS_solution(temp_sol, U_sol[itime], P_sol[itime]);
    }
}

void PHYSICS::eval_mean_solution_over_time(MATRIX &solution_over_times, VECTOR &time_avg_sol)
{
    time_avg_sol.complete_reset();
    time_avg_sol.initialize(solution_over_times.nCol);
    if (isStationary)
    {
        time_avg_sol = solution_over_times.get_row(0);
    }
    else
    {
        for (int itime = 0; itime < solution_times.length-1; itime++)
        {
            VECTOR time_sol_i = solution_over_times.get_row(itime);
            VECTOR time_sol_j = solution_over_times.get_row(itime+1);
            VECTOR mean_time_sol = (time_sol_i+time_sol_j) / 2.0;
            time_avg_sol += mean_time_sol * solution_deltaT[itime];
        }
        time_avg_sol /= solution_times.get_last() - solution_times[0];
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
    element_value.complete_reset();
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
    element_value.complete_reset();
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
    centroids_value.complete_reset();
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
    smoothed_value.complete_reset();
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
    smoothed_value.complete_reset();
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

// invert std::vector<VECTOR> sizes
void PHYSICS::invert_sizes(std::vector<VECTOR> &vec, std::vector<VECTOR> &inverted_vec)
{
    inverted_vec.clear();
    int n_comp = int(vec.size());
    int n_val = vec[0].length;
    inverted_vec.resize(n_val);
    for (int ival = 0; ival < n_val; ival++)
    {
        inverted_vec[ival].setZeros(n_comp);
    }
    for (int icomp = 0; icomp < n_comp; icomp++)
    {
        for (int ival = 0; ival < n_val; ival++)
        {
            inverted_vec[ival][icomp] = vec[icomp][ival];
        }
    }
}

void PHYSICS::eval_gradient(VECTOR &value, std::vector<VECTOR> &gradient, int eval_method)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[nodal_value]
    // gradient structure: gradient[node][gradient_component]
    //eval_method = 0: Lampeato. Other cases not implemented
    int n_val = value.length;
    gradient.clear();
    gradient.resize(n_val);
    if (eval_method == 0) // Lampeato
    {
        VECTOR weights; 
        weights.setZeros(nNodes_v);
        for (int inod = 0; inod < n_val; inod++)
        {
            gradient[inod].complete_reset();
            gradient[inod].setZeros(dim);
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
                    gradient[iglob][icomp] += element_gradient_value[icomp] * temp_weight;
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
        for (int inod = 0; inod < n_val; inod++)
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                // if ((iglob == 2392))
                // {
                //     std::cout << "\npre-comp " << icomp << ": " << gradient[icomp][iglob] << "\n";
                // }
                gradient[inod][icomp] /= weights[inod];
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

void PHYSICS::eval_gradient(MATRIX &value, std::vector<MATRIX> &gradient, int eval_method)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[component][nodal_value]
    // gradient structure: gradient[node][value_component][derivative_component]
    //initialize solution
    int n_val = value.nCol;
    int n_comps = value.nRow;
    gradient.clear();
    gradient.resize(n_val);
    for (int inod = 0; inod < n_val; inod++)
    {
        gradient[inod].complete_reset();
    }
    
    for (int icomp = 0; icomp < n_comps; icomp++)
    {
        VECTOR comp_value = value.get_row(icomp);
        std::vector<VECTOR> temp_grad;
        eval_gradient(comp_value, temp_grad, eval_method);
        for (int inod = 0; inod < n_val; inod++)
        {
            VECTOR temp_nodal_grad = temp_grad[inod];
            gradient[inod].append_row(temp_nodal_grad);
        }
    }
}

void PHYSICS::eval_divergence(VECTOR &value, VECTOR &divergence)
{
    int n_val = value.length;
    std::vector<VECTOR> gradient;
    divergence.complete_reset();
    divergence.setZeros(n_val);
    eval_gradient(value, gradient);
    for (int inod = 0; inod < n_val; inod++)
    {
        divergence[inod] = gradient[inod].sum_comps();
    }
}

void PHYSICS::eval_divergence(MATRIX &value, std::vector<VECTOR> &divergence)
{
    int n_val = value.nCol;
    int n_comps = value.nRow;
    divergence.resize(n_val);
    for (int inod = 0; inod < n_val; inod++)
    {
        divergence[inod].complete_reset();
        divergence[inod].setZeros(n_comps);
    }
    std::vector<MATRIX> gradient;
    eval_gradient(value, gradient);
    for (int inod = 0; inod < n_val; inod++)
    {
        MATRIX temp_grad = gradient[inod];
        for (int icomp = 0; icomp < n_comps; icomp++)
        {
            VECTOR temp_comp_grad = temp_grad.get_row(icomp);
            divergence[inod][icomp] = temp_comp_grad.sum_comps();
        }
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
    dir_gradient.complete_reset();
    dir_gradient.setZeros(nodes.length);
    std::vector<VECTOR> gradient;
    eval_gradient(value, gradient);
    for (int inod = 0; inod < nodes.length; inod++)
    {
        int iglob = nodes[inod];
        VECTOR temp_dir = directions[inod];
        dir_gradient[inod] = gradient[iglob].dot(directions[inod]);
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
    // dir_gradient structure: gradient[node][value_component_derivative]
    // N.B.: specific node means it belongs to [0,nodes-length-1]
    int n_comps = value.nRow;
    dir_gradient.resize(nodes.length);
    for (int inod = 0; inod < nodes.length; inod++)
    {
        dir_gradient[inod].complete_reset();
        dir_gradient[inod].setZeros(n_comps);
    }
    for (int icomp = 0; icomp < n_comps; icomp++)
    {
        VECTOR comp_value = value.get_row(icomp);
        VECTOR temp_dir_grad;
        eval_directional_gradient(comp_value, nodes, directions, temp_dir_grad);
        for (int inod = 0; inod < nodes.length; inod++)
        {
            dir_gradient[inod][icomp] = temp_dir_grad[inod];
        }
    }
}

void PHYSICS::eval_gradient_norm(std::vector<VECTOR> &gradient, VECTOR &norm)
{
    norm.complete_reset();
    norm.setZeros(nNodes_v);
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        for (int icomp = 0; icomp < dim; icomp++)
        {
            norm[inod] += gradient[inod][icomp]*gradient[inod][icomp];
        }
    }
    norm.squared_root();
}

void PHYSICS::eval_gradient_norm(std::vector<MATRIX> &gradient, VECTOR &norm)
{
    int n_val = int(gradient.size());
    norm.complete_reset();
    norm.setZeros(n_val);
    for (int inod = 0; inod < n_val; inod++)
    {
        MATRIX nodal_gradient = gradient[inod];
        norm[inod] = nodal_gradient.normFro();
    }
}

void PHYSICS::eval_WSS(MATRIX &value, VECTOR_INT &nodes, std::vector<VECTOR> &normals, std::vector<VECTOR> &WSS)
{
    // suppose the value interpolated in the coord_v
    // value structure: value[component][nodal_value]
    // normals structure: normals[node][normal_component]
    // WSS structure: 
    //                  WSS[0][node] WSS at node
    //                  WSS[i!=0][node] mu*(grad(U)*normal)*tangent_i ("ith tangential component contributing to WSS")
    int n_nodes = nodes.length;
    if (n_nodes != int(normals.size()))
    {
        throw_line("ERROR: the nodes and the normals for the WSS evaluation are defined in a different number of points, or in a uncompatible format\n");
    }

    std::vector<VECTOR> normal_gradient;
    eval_directional_gradient(value, nodes, normals, normal_gradient);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        bool is_nan = normal_gradient[icomp].check_nan();
        if (is_nan)
        {
            std::cout << "comp: " << icomp << "\n";
            throw_line("ERROR: nan value eval_directional_gradient for the WSS\n");
        }
    }
    int n_tan = dim-1;
    std::vector<std::vector<VECTOR>> tangent_vectors;
    // build tangent vectors
    build_tangents_from_normals(normals, tangent_vectors);

    WSS.clear();
    WSS.resize(dim);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        WSS[icomp].complete_reset();
        WSS[icomp].setZeros(n_nodes);
    }
    // project normal gradient into tangent direction(2D)/plane(3D)
    for (int inod = 0; inod < n_nodes; inod++)
    {
        VECTOR nodal_grad = normal_gradient[inod]; 
        for (int itan = 0; itan < n_tan; itan++)
        {
            prec temp_grad_proj = mu * nodal_grad.dot(tangent_vectors[inod][itan]); // containts the mu for the WSS
            WSS[itan+1][inod] = temp_grad_proj*temp_grad_proj;
        }        
    }
    for (int itan = 0; itan < n_tan; itan++)
    {
        WSS[0] += WSS[itan+1];
    }
    for (int icomp = 0; icomp < dim; icomp++)
    {
        WSS[icomp].squared_root();
    }
}

void PHYSICS::eval_WSS(MATRIX &value, int bound_id, VECTOR normal, std::vector<VECTOR> &WSS)
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

void PHYSICS::eval_WSS_avg_over_time(std::vector<MATRIX> &value, VECTOR_INT &nodes, std::vector<VECTOR> &normals, std::vector<VECTOR> &WSS)
{
    // WSS structure: 
    //                  WSS[0][node] WSS at node
    //                  WSS[i!=0][node] (grad(U)*normal)*tangent_i ("ith tangential component contributing to WSS")
    WSS.clear();
    WSS.resize(dim);
    std::vector<MATRIX> WSS_over_times(dim);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        WSS_over_times[icomp].complete_reset();
    }
    for (int itime = 0; itime < int(value.size()); itime++)
    {
        bool is_nan = value[itime].check_nan(true);
        std::vector<VECTOR> temp_WSS;
        eval_WSS(value[itime], nodes, normals, temp_WSS);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            VECTOR WSS_comp = temp_WSS[icomp];
            is_nan = WSS_comp.check_nan();
            if (is_nan)
            {
                std::cout << "\ntime: " << itime << " contains nan\n\n";
                throw_line("ERROR in eval_WSS_avg_over_time\n");
            }
            WSS_over_times[icomp].append_row(WSS_comp);
        }
        
    }

    // integral average over times
    for (int icomp = 0; icomp < dim; icomp++)
    {
        eval_mean_solution_over_time(WSS_over_times[icomp], WSS[icomp]);
    } 
}

void PHYSICS::build_tangents_from_normal(VECTOR &normal, std::vector<VECTOR> &tangent_vectors)
{
    // normal structure: normal[component]
    // tangents structure: tangent_vectors[tangent_id][tangent_component]

    if (normal.length != dim)
    {
        throw_line("ERROR: normal has a size different from the dimension\n");
    }
    tangent_vectors.clear();
    int n_tan = dim-1;
    tangent_vectors.resize(n_tan);
    for (int icomp = 0; icomp < n_tan; icomp++)
    {
        tangent_vectors[icomp].complete_reset();
        tangent_vectors[icomp].setZeros(dim);
    }
    if (normal.norm() > 1e-12)
    {
        switch (n_tan)
        {
            case 1: // 2D
            {
                // tangent_vectors[inod][0].setZeros(dim);

                tangent_vectors[0][0] = -normal[1];
                tangent_vectors[0][1] = normal[0];
                // normalize the tangent vector
                tangent_vectors[0] /= tangent_vectors[0].norm();
                break;
            }
            case 2: // 3D
            {
                // first tangent vector
                tangent_vectors[0][0] = -normal[1];
                tangent_vectors[0][1] = normal[0];
                // normalize the tangent vector
                prec first_tan_norm = tangent_vectors[0].norm();
                tangent_vectors[0] /= first_tan_norm;

                // second tangent vector, orthogonal to the first one
                tangent_vectors[1][0] = -normal[2] + (normal[1]*normal[1]*normal[2])/(first_tan_norm*first_tan_norm);
                tangent_vectors[1][1] = -(normal[0]*normal[1]*normal[2])/(first_tan_norm*first_tan_norm);
                tangent_vectors[1][2] = normal[0];
                tangent_vectors[1] /= tangent_vectors[1].norm();
                break;
            }   
            default:
            {
                throw_line("ERROR with dimensions and number of tangents\n");
                break;
            }
        }
    }
}
    
void PHYSICS::build_tangents_from_normals(std::vector<VECTOR> &normals, std::vector<std::vector<VECTOR>> &tangent_vectors)
{
    // normals structure: normals[node][component]
    // tangents structure: tangent_vectors[node][tangent_id][tangent_component]

    int n_nodes = int(normals.size());
    tangent_vectors.clear();
    tangent_vectors.resize(n_nodes);
    for (int inod = 0; inod < n_nodes; inod++)
    {
        VECTOR nodal_normal = normals[inod];
        build_tangents_from_normal(nodal_normal, tangent_vectors[inod]); 
    }
}

void PHYSICS::print_mesh_for_mmg(std::string file_name, std::string rel_path)
{
    if (dim != 2) throw_line("NOT IMPLEMENTED REMESHING FOR THE 3D CASE");

    // creat mesh file
    std::string currFilePath;
    if (rel_path == "")
    {
        currFilePath = file_name;
    }
    else
    {
        currFilePath = rel_path + "/" + file_name;
    }
    std::ofstream file;
    file.open(currFilePath);

    // initialize mesh
    file << "MeshVersionFormatted 2\n";
    file << "\n";
    file << "\n";
    file << "Dimension " << dim << "\n";
    file << "\n";
    file << "\n";

    // write nodes
    file << "Vertices\n";
    file << nNodes << "\n";
    for (int inod = 0; inod < nNodes; inod++)
    {
        prec* temp_coord = coord[inod];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            prec temp_val = temp_coord[icomp];
            if (abs(temp_val) < 1e-14) temp_val = 0.0;
            file << temp_val << " ";
        }
        if (dim == 2)
        {
            file << "0 ";
        }
        file << "\n";
    }
    
    file << "\n";
    file << "\n";

    // write elements
    file << "Triangles\n";
    file << nElem << "\n";
    for (int iel = 0; iel < nElem; iel++)
    {
        int* temp_nodes = elem[iel];
        for (int inod = 0; inod < dim+1; inod++)
        {
            file << (temp_nodes[inod]+1) << " ";
        }
        file << (elem_geo_entities_ids_v[4*iel]);
        file << "\n";
    }

    file << "\n";
    file << "\n";
    file << "End\n";

    // close mesh file
    file.close();
}
//--------------------------------------------
// END PHYSICS CLASS
//--------------------------------------------

