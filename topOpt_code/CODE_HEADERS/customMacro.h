#pragma once

#include "codeHeader.h"

void freeVecOfVECTOR(std::vector<VECTOR> vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        // vec[i].dlt();
    }
    std::vector<VECTOR>().swap(vec);
}

