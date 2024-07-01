#include "scrape_meshes.h"

int main(int argc, char *argv[])
{
    // bool BINREAD;
    if (argc < 2) throw_line("ERROR: missing file input name");
    // if (argc < 3) BINREAD = false;
    // else BINREAD = argv[2];
    char* file_name = argv[1];
    std::string input_file = file_name;
    MATRIX coord;
    MATRIX_INT elem;
    VECTOR_INT elem_id;
    int dim = 2;
    MESH_SCRAPER scraper;

    scraper.extract_mesh_from_msh(input_file, dim, coord, elem, elem_id);

    std::string output_text_mph = scraper.print_mesh_in_mphtxt(coord, elem, elem_id, dim);
    // Create and open a text file
    std::string output_file_mph = "output.mphtxt";
    std::ofstream output_stream_mph(output_file_mph);
    output_stream_mph << output_text_mph;
    output_stream_mph.close();

    std::string output_text_msh = scraper.print_mesh_in_msh(coord, elem, elem_id, dim);
    // Create and open a text file
    std::string output_file_msh = "output.msh";
    std::ofstream output_stream_msh(output_file_msh);
    output_stream_msh << output_text_msh;
    output_stream_msh.close();

    std::string output_text_mesh = scraper.print_mesh_in_mesh(coord, elem, elem_id, dim);
    // Create and open a text file
    std::string output_file_mesh = "output.mesh";
    std::ofstream output_stream_mesh(output_file_mesh);
    output_stream_mesh << output_text_mesh;
    output_stream_mesh.close();
    
    return 0;
}