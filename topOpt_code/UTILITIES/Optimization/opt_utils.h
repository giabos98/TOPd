#include "../../CODE_HEADERS/codeHeader.h"

class TIME_PROFILER
{
    public:

    //-----------------
    // PROPERTIES
    //-----------------
    VECTOR times;

    //-----------------
    // METHODS
    //-----------------
    TIME_PROFILER() 
    {
        prec time = omp_get_wtime();
        times.append(time);
    }

    void initialize();
    void set();
    void profile(std::string time_name);
};