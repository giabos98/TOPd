#include "opt_utils.h"

void TIME_PROFILER::initialize()
{
    set();
}

void TIME_PROFILER::set()
{
    prec time = omp_get_wtime();
    times.append(time);
}

void TIME_PROFILER::profile(std::string time_name)
{
    set();
    std::cout << "\n--| " << time_name << ": " << (times.get(-1) - times.get(-2)) << "\n";
}