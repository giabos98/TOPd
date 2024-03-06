#include "../CODE_HEADERS/codeHeader.h"

std::string scrape_3D_mesh_into_2D(std::string input_file)
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


int main(int argc, char *argv[])
{
    bool BINREAD;
    if (argc < 2) throw std::runtime_error("usage: missing file input name");
    if (argc < 3) BINREAD = false;
    else BINREAD = argv[2];
    char* file_name = argv[1];
    std::string input_file = file_name;
    std::string res_text = scrape_3D_mesh_into_2D(input_file);
    // Create and open a text file
    std::string output_file = "output.mphtxt";
    std::ofstream output_stream(output_file);
    // Write to the file
    output_stream << res_text;
    // Close the file
    output_stream.close();
    return 0;
}