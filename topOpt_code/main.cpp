#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"
#include "TopOpt.h"
int PARALLEL::nThread = 4;
int main()//int argc, char* argv[])
{    
    // char* n_threads_char = argv[1];
    // char* n_reps_char = argv[2];
    // std::string n_threads_str = n_threads_char;
    // std::string n_reps_str = n_reps_char;
    // int n_threads = std::stoi(n_threads_str);
    // int n_reps = std::stoi(n_reps_str);
    VECTOR_INT n_th(1);
    n_th[0] = 20;
    // n_th[1] = 2;
    // n_th[2] = 4;
    // n_th[3] = 8;
    // n_th[4] = 14;
    // n_th[5] = 20;
    int n_reps = 5;

    for (int ith = 0; ith < n_th.length; ith++)
    {
        int n_threads = n_th[ith];

        prec startTime = omp_get_wtime();
        VECTOR general_times(4);

        std::string inputFileNS = "./INPUT_FILES/readProblemNS.txt";
        TOP_OPT topOpt(inputFileNS, n_threads, n_reps, general_times);

        //*-*-*--*-*-*-*-*-*-*
        // topOpt.solve();
        //*-*-*-*-*-*-*-*-*-*-*-

        // prec endTime = omp_get_wtime();
        // prec totalTime = endTime-startTime;

        // topOpt.print_stats(totalTime);
        general_times[0] = general_times[0] - startTime;
        general_times[3] = omp_get_wtime() - startTime;

        std::cout << "\n-| total time: " << general_times[3] << "\n";
        std::string times_file = "./results/porosity_test/Capri_HPC_test/general_times_" + std::to_string(n_threads) + ".txt";
        general_times.printFile(times_file.c_str());
    }

    std::cout << "\n---| END |---\n";

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