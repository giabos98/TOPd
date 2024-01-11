#include "../CODE_HEADERS/codeHeader.h"
#include "loadFromComsolFile.h"

int main(int argc, char *argv[])
{
    bool BINREAD;
    if (argc < 2) throw std::runtime_error("usage: missing file input name");
    if (argc < 3) BINREAD = false;
    else BINREAD = argv[2];

    char* fileName = argv[1];
    int N = strlen(fileName);

    const char* folderPath = "MESHES_COMSOL/";

    loadMesh(fileName, N, folderPath, BINREAD);
};
