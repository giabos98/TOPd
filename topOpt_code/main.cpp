#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"
#include "TopOpt.h"
int PARALLEL::nThread = 4;
int main()
{    
    prec startTime = omp_get_wtime();

    std::string inputFileNS = "INPUT_FILES/readProblemNS.txt";
    TOP_OPT topOpt(inputFileNS);
    // topOpt.physics.print_mesh_for_mmg("prova.mesh");

    //*-*-*--*-*-*-*-*-*-*
    topOpt.solve();
    //*-*-*-*-*-*-*-*-*-*-*-

    prec endTime = omp_get_wtime();

    prec totalTime = endTime - startTime;
    topOpt.print_stats(totalTime);

    return 0;

}