#include "scrape_meshes.h"

void MESH_SCRAPER::initialize()
{
    // empty constructor
}

std::string MESH_SCRAPER::scrape_3D_mesh_into_2D(std::string input_file) // BETA
{
    std::cout << "\n" << "---| SCRAPE 3D COMSOL MESH INTO 2D |---" << "\n\n";
    std::string output_text = "";
    std::ifstream mesh_stream;
    mesh_stream.open(input_file);
    std::string line;
    std::istringstream iss;

    for (int i = 0; i < 6; i++)
    {
        getline(mesh_stream, line);
        output_text += line + "\n";
    }
    STREAM::getLines(mesh_stream, line, 1);
    output_text += "5 mesh2 \n";
    for (int i = 0; i < 9; i++)
    {
        getline(mesh_stream, line);
        output_text += line + "\n";
    }
    STREAM::getLines(mesh_stream, line, 1);
    output_text += "2 # sdim \n";
    int n_nodes = 0;
    STREAM::getValue(mesh_stream, line, iss, n_nodes);
    output_text += std::to_string(n_nodes) + " # number of mesh points \n";
    STREAM::getLines(mesh_stream, line, 1);
    output_text += "0 # lowest mesh point index \n";
    getline(mesh_stream, line);
    output_text += line + "\n";
    getline(mesh_stream, line);
    output_text += "# Mesh point coordinates \n";
    for (int i = 0; i < n_nodes; i++)
    {
        VECTOR temp_coord(2);
        STREAM::getRowVector(mesh_stream, line, iss, temp_coord);
        output_text += std::to_string(temp_coord[0]) + " " + std::to_string(temp_coord[1]) + "\n";
    }
    for (int i = 0; i < 9; i++)
    {
        getline(mesh_stream, line);
        output_text += line + "\n";
    }
    int n_b_el = 0;
    STREAM::getValue(mesh_stream, line, iss, n_b_el);
    output_text += std::to_string(n_b_el) + " # number of elements \n";
    getline(mesh_stream, line);
    output_text += line + "\n";
    for (int i = 0; i < n_b_el; i++)
    {
        VECTOR_INT temp_coord(2);
        STREAM::getRowVector(mesh_stream, line, iss, temp_coord);
        output_text += std::to_string(temp_coord[0]) + " " + std::to_string(temp_coord[1]) + "\n";
    }
    for (int i = 0; i < 9; i++)
    {
        getline(mesh_stream, line);
        output_text += line + "\n";
    }
    int n_el = 0;
    STREAM::getValue(mesh_stream, line, iss, n_el);
    output_text += std::to_string(n_el) + " # number of elements \n";
    getline(mesh_stream, line);
    output_text += line + "\n";
    for (int i = 0; i < n_el; i++)
    {
        VECTOR_INT temp_coord(3);
        STREAM::getRowVector(mesh_stream, line, iss, temp_coord);
        output_text += std::to_string(temp_coord[0]) + " " + std::to_string(temp_coord[1]) + " " + std::to_string(temp_coord[2]) +"\n";
    }
    for (int i = 0; i < 2; i++)
    {
        getline(mesh_stream, line);
        output_text += line + "\n";
    }
    mesh_stream.close();
    return output_text;
}

//-----------------------------------------------------------------------------------------
// EXTRACT COORD, ELEM AND ELEM_ID FROM A .MSH MESH FILE (version 2.2)
//-----------------------------------------------------------------------------------------
void MESH_SCRAPER::extract_mesh_from_msh(std::string input_file, int dim, MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id)
{
    std::cout << "\n" << "---| EXTRACT COORD, ELEM AND ELEM_ID FROM A .MPH MESH FILE |---" << "\n\n";
    std::ifstream mesh_stream;
    mesh_stream.open(input_file);
    std::string line;
    std::istringstream iss;
    STREAM::getLines(mesh_stream, line, 4);
    int n_nodes;
    STREAM::getValue(mesh_stream, line, iss, n_nodes);
    coord.complete_reset();
    coord.initialize(n_nodes, dim);
    int n_info;
    n_info = dim+1;
    for (int inod = 0; inod < n_nodes; inod++)
    {
        VECTOR temp_val(n_info);
        STREAM::getRowVector(mesh_stream, line, iss, temp_val);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            coord[inod][icomp] = temp_val[icomp+1];
        }
    }
    STREAM::getLines(mesh_stream, line, 2);
    int n_max_el;
    STREAM::getValue(mesh_stream, line, iss, n_max_el);
    elem.complete_reset();
    int n_nod_x_el = dim+1;
    elem.initialize(n_max_el, n_nod_x_el);
    elem_id.complete_reset();
    elem_id.initialize(n_max_el);
    int not_node_info = 5;
    n_info = not_node_info + n_nod_x_el;
    int count_el = 0;
    for (int iel = 0; iel < n_max_el; iel++)
    {
        VECTOR_INT temp_val(n_info);
        STREAM::getRowVector(mesh_stream, line, iss, temp_val);
        if (int(temp_val[1]) == 2)
        {
            elem_id[count_el] = temp_val[not_node_info-1]-1;
            for (int inod = 0; inod < n_nod_x_el; inod++)
            {
                elem[count_el][inod] = temp_val[not_node_info+inod]-1;
            }
            count_el++;
        }
    }
    elem.shrinkRows(count_el);
    elem_id.shrink(count_el);
}

//-----------------------------------------------------------------------------------------
// CREATE THE STRING TO BE PRINTED IN A MPHTXT FILE TO READ IT IN COMSOL
//-----------------------------------------------------------------------------------------
std::string MESH_SCRAPER::print_mesh_in_mphtxt(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim)
{
    std::cout << "\n" << "---| BUILD .MPHTXT MESH FROM COORDS AND ELEMS (with elem domain id) |---" << "\n\n";
    std::string text = "";

    text += "# Major & minor version\n";
    text += "0 1\n";
    text += "1 # number of tags\n";
    text += "# Tags\n"; 
    text += "5 mesh1\n "; 
    text += "1 # number of types\n"; 
    text += "# Types\n"; 
    text += "3 obj\n "; 
    text += "\n";
    text += "# --------- Object 0 ----------\n"; 
    text += "\n"; 
    text += "0 0 1\n"; 
    text += "4 Mesh # class\n"; 
    text += "4 # version\n";
    text += std::to_string(dim) + " # sdim\n";
    int n_nodes = coord.nRow;
    text += std::to_string(n_nodes) + " # number of mesh vertices\n";
    text += "0 # lowest mesh vertex index\n";
    text += "\n";
    text += "# Mesh vertex coordinates\n";
    for (int inod = 0; inod < n_nodes; inod++)
    {
        std::string temp_line = "";
        prec* temp_coord = coord[inod];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            temp_line += std::to_string(temp_coord[icomp]) + " ";
        }
        text += temp_line + "\n";
    }
    text += "\n";
    text += "1 # number of element types\n";
    text += "\n";
    text += "# Type #0";
    text += "\n";
    text  += "3 tri # type name";
    text += "\n";
    text += "\n";
    int n_nod_x_el = dim+1;
    text += std::to_string(n_nod_x_el) + " # number of vertices per element\n";
    int n_el = elem.nRow;
    text += std::to_string(n_el) + " # number of elements\n";
    text += "# Elements\n";
    for (int iel = 0; iel < n_el; iel++)
    {
        std::string temp_line = "";
        int* temp_nodes = elem[iel];
        for (int inod = 0; inod < n_nod_x_el; inod++)
        {
            temp_line += std::to_string(temp_nodes[inod]) + " ";
        }
        text += temp_line + "\n";
    }
    text += "\n";
    int n_el_id = elem_id.length;
    if (n_el_id != n_el) throw_line("ERROR: incompatible number of elements and domain id vector length\n");
    text += std::to_string(n_el_id) + " # number of geometric entity indices\n";
    text += "# Geometric entity indices\n";
    for (int iel = 0; iel < n_el_id; iel++)
    {
        text += std::to_string(elem_id[iel]) + "\n";
    } 
    return text;
}

//-----------------------------------------------------------------------------------------
// CREATE THE STRING TO BE PRINTED IN A .MSH FILE TO READ IT IN GMSH
//-----------------------------------------------------------------------------------------
std::string MESH_SCRAPER::print_mesh_in_msh(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim)
{
    std::cout << "\n" << "---| BUILD .MSH MESH FROM COORDS AND ELEMS (with elem domain id) |---" << "\n\n";
    std::string text = "";

    text += "$MeshFormat\n";
    text += "2.2 0 8\n";
    text += "$EndMeshFormat\n";
    text += "$Nodes\n";
    int n_nodes = coord.nRow;
    text += std::to_string(n_nodes) + "\n"; 
    for (int inod = 0; inod < n_nodes; inod++)
    {
        std::string temp_text = std::to_string(inod+1) + " ";
        prec* temp_coord = coord[inod];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            temp_text += std::to_string(temp_coord[icomp]) + " ";
        }
        if (dim == 2)
        {
            temp_text += "0 ";
        }
        text += temp_text + "\n";
    }
    text += "$EndNodes\n";
    text += "$Elements\n";
    int n_elem = elem.nRow;
    int n_nodes_x_el = elem.nCol;
    text += std::to_string(n_elem) + "\n";
    for (int iel = 0; iel < n_elem; iel++)
    {
        std::string temp_text = std::to_string(iel+1) + " ";
        temp_text += "2 2 0 ";
        temp_text += std::to_string(elem_id[iel]+1) + " ";
        int* temp_nodes = elem[iel];
        for (int inod = 0; inod < n_nodes_x_el; inod++)
        {
            temp_text += std::to_string(temp_nodes[inod]+1) + " ";
        }
        text += temp_text + "\n";
    }
    text += "$EndElements\n";
    return text;
}

//-----------------------------------------------------------------------------------------
// CREATE THE STRING TO BE PRINTED IN A .MESH FILE TO READ IT IN GMSH
//-----------------------------------------------------------------------------------------
std::string MESH_SCRAPER::print_mesh_in_mesh(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim)
{
    std::cout << "\n" << "---| BUILD .MESH MESH FROM COORDS AND ELEMS (with elem domain id) |---" << "\n\n";
    std::string text = "";

    text += "MeshVersionFormatted 2\n";
    text += "\n";
    text += "\n";
    text += "Dimension " + std::to_string(dim) + "\n";
    text += "\n";
    text += "\n";
    text += "Vertices\n";
    int n_nodes = coord.nRow;
    text += std::to_string(n_nodes) + "\n"; 
    for (int inod = 0; inod < n_nodes; inod++)
    {
        std::string temp_text = "";
        prec* temp_coord = coord[inod];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            temp_text += std::to_string(temp_coord[icomp]) + " ";
        }
        if (dim == 2)
        {
            temp_text += "0 ";
        }
        text += temp_text + "\n";
    }
    text += "\n";
    text += "\n";
    text += "Triangles\n";
    int n_elem = elem.nRow;
    int n_nodes_x_el = elem.nCol;
    text += std::to_string(n_elem) + "\n";
    for (int iel = 0; iel < n_elem; iel++)
    {
        std::string temp_text = "";
        int* temp_nodes = elem[iel];
        for (int inod = 0; inod < n_nodes_x_el; inod++)
        {
            temp_text += std::to_string(temp_nodes[inod]+1) + " "; // indices starting from 1
        }
        temp_text += std::to_string(elem_id[iel]+1);
        text += temp_text + "\n";
    }
    text += "\n";
    text += "\n";
    text += "End\n";
    return text;
}

//-----------------------------------------------------------------------------------------
// CREATE THE STRING TO PRINT A SCALAR SOLUTION IN A .SOL FILE
//-----------------------------------------------------------------------------------------
std::string MESH_SCRAPER::print_scalar_solution_in_sol(VECTOR &solution, int dim)
{
    int nNodes = solution.length;
    std::string sol_string = "";
    sol_string += "MeshVersionFormatted 2\n";
    sol_string += "\n";
    sol_string += "Dimension " + std::to_string(dim) + "\n";
    sol_string += "\n";
    sol_string += "SolAtVertices\n";
    sol_string += std::to_string(nNodes) + "\n";
    sol_string += "1 1\n";
    sol_string += "\n";
    for (int inod = 0; inod < nNodes; inod++)
    {
        sol_string += std::to_string(solution[inod]) + "\n";
    }
    sol_string += "\n";
    sol_string += "End\n";

    return sol_string;
}