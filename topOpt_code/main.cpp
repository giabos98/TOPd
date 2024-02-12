#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"
#include "TopOpt.h"
int PARALLEL::nThread = 4;
int main()
{    
    prec startTime = omp_get_wtime();

    std::string inputFileNS = "./INPUT_FILES/readProblemNS.txt";
    TOP_OPT topOpt(inputFileNS);

    //*-*-*--*-*-*-*-*-*-*
    topOpt.solve();
    //*-*-*-*-*-*-*-*-*-*-*-

    prec endTime = omp_get_wtime();

    prec totalTime = endTime - startTime;
    topOpt.print_stats(totalTime);

    // int nth = 20;
    // int n = 1e8;
    // prec startTime; prec endTime;
    // prec res;

    // res = 0;
    // VECTOR vec1(n);
    // VECTOR vec2(n);
    // for (int i = 0; i < n; i++)
    // {
    //     vec1[i] = 1;
    //     vec2[i] = 1;
    // }
    // VECTOR vec3(n);

    // omp_set_num_threads(nth);
    // startTime = omp_get_wtime();
    // // #pragma omp parallel for
    // // for (int i = 0; i < n; i++)
    // // {
    // //     res += vec1[i] * vec2[i];
    // // }
    // vec1.sum(vec2.P, vec3.P);

    // endTime = omp_get_wtime();

    // std::cout << "\ncstm: " << endTime-startTime << "\n";

    // res = 0;
    // std::vector<double> v1(n);
    // std::vector<double> v2(n);
    // for (int i = 0; i < n; i++)
    // {
    //     v1[i] = 1;
    //     v2[i] = 1;
    // }
    // std::vector<prec> v3(n);

    // omp_set_num_threads(nth);
    // startTime = omp_get_wtime();
    // // #pragma omp parallel for
    // for (int i = 0; i < n; i++)
    // {
    //     v3[i] = v1[i] + v2[i];
    // }
    // endTime = omp_get_wtime();
    // std::cout << "\nstd: " << endTime-startTime << "\n";

    // res = 0;
    // prec* p1 = (prec*) malloc(n*sizeof(prec));
    // prec* p2 = (prec*) malloc(n*sizeof(prec));
    // for (int i = 0; i < n; i++)
    // {
    //     p1[i] = 1;
    //     p2[i] = 1;
    // }
    // prec* p3 = (prec*) malloc(n*sizeof(prec));

    // omp_set_num_threads(nth);
    // startTime = omp_get_wtime();
    // // #pragma omp parallel for
    // for (int i = 0; i < n; i++)
    // {
    //     p3[i] = p1[i] + p2[i];
    // }
    // endTime = omp_get_wtime();
    // std::cout << "\nptr: " << endTime-startTime << "\n";

    return 0;

}