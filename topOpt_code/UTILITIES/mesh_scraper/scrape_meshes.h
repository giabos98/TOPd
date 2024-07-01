#include "../../CODE_HEADERS/codeHeader.h"

class MESH_SCRAPER
{
    public:

    //-----------------
    // PROPERTIES
    //-----------------

    //-----------------
    // METHODS
    //-----------------
    MESH_SCRAPER() {}

    void initialize();

    //-----------------------------------------------------------------------------------------
    // SCRAPE 3D COMSOL MESHES INTO 2D ONES
    //-----------------------------------------------------------------------------------------
    std::string scrape_3D_mesh_into_2D(std::string input_file); // BETA

    //-----------------------------------------------------------------------------------------
    // EXTRACT COORD, ELEM AND ELEM_ID FROM A .MSH MESH FILE (version 2.2)
    //-----------------------------------------------------------------------------------------
    void extract_mesh_from_msh(std::string input_file, int dim, MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id);

    //-----------------------------------------------------------------------------------------
    // CREATE THE STRING TO BE PRINTED IN A MPHTXT FILE TO READ IT IN COMSOL
    //-----------------------------------------------------------------------------------------
    std::string print_mesh_in_mphtxt(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim);

    //-----------------------------------------------------------------------------------------
    // CREATE THE STRING TO BE PRINTED IN A .MSH FILE TO READ IT IN GMSH
    //-----------------------------------------------------------------------------------------
    std::string print_mesh_in_msh(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim);

    //-----------------------------------------------------------------------------------------
    // CREATE THE STRING TO BE PRINTED IN A .MESH FILE TO READ IT IN GMSH
    //-----------------------------------------------------------------------------------------
    std::string print_mesh_in_mesh(MATRIX &coord, MATRIX_INT &elem, VECTOR_INT &elem_id, int dim);

    //-----------------------------------------------------------------------------------------
    // CREATE THE STRING TO PRINT A SCALAR SOLUTION IN A .SOL FILE
    //-----------------------------------------------------------------------------------------
    std::string print_scalar_solution_in_sol(VECTOR &solution, int dim);

};