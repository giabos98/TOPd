#include "CODE_HEADERS/codeHeader.h"
#include "darcyProblem.h"

void PROBLEM_DARCY::localBasis()
{
    std::cout << "\n----------\n--| LOCAL BASIS  |--\n----------\n";
    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    std::shared_ptr<prec*[]>  coord = (*physics).coord.PP; 
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    (*physics).Volume.initialize(nElem);
    std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
    // int nElem_v = (*physics).nElem_v;
    // std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
    // std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    // (*physics).Volume_v.initialize(nElem_v);
    // std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    
    (*physics).Bloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Bloc = (*physics).Bloc.PP;
    (*physics).Cloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Cloc = (*physics).Cloc.PP;

    // VELOCITY
    // (*physics).Bloc_v.initialize(nElem_v, dim+1);
    // std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    // (*physics).Cloc_v.initialize(nElem_v, dim+1);
    // std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    

    switch (dim)
    {
        case 2: // 3 nodes for triangles
        {
            // PRESSURE NODES
            for (int iel = 0; iel < nElem; iel++)
            {
                int* locNodes = elem[iel];
                prec p0x = coord[locNodes[0]][0]; prec p0y = coord[locNodes[0]][1];
                prec p1x = coord[locNodes[1]][0]; prec p1y = coord[locNodes[1]][1];
                prec p2x = coord[locNodes[2]][0]; prec p2y = coord[locNodes[2]][1];

                prec detVand = p2y*(p1x-p0x) - p1y*(p2x-p0x) + p0y*(p2x-p1x);
                Volume[iel] = std::abs(detVand)/2;
                for (int inod = 0; inod < 3; inod++)
                {
                    int n1 = (inod+1)%3;
                    int n2 = (inod+2)%3;
                    Bloc[iel][inod] = (coord[locNodes[n1]][1] - coord[locNodes[n2]][1])/detVand;
                    Cloc[iel][inod] = (coord[locNodes[n2]][0] - coord[locNodes[n1]][0])/detVand;
                }
            }
            break;
        }
        case 3: // 4 nodes for tetrahedrums
        {
            // PRESSURE
            (*physics).Dloc.initialize(nElem, 4);
            std::shared_ptr<prec*[]> Dloc = (*physics).Dloc.PP;
            for (int iel = 0; iel < nElem; iel++)
            {
                int* locNodes = elem[iel];
                MATRIX Vand(4,dim);
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        Vand[i][j] = coord[locNodes[i]][j];
                    }
                }

                prec detVand = 0;
                //---
                prec temp_a; 
                prec temp_b;
                prec temp_c;
                prec temp_d;
                int sign = 1;
                for (int inod = 0; inod < 4; inod++)
                {
                    int n1 = (inod+1)%4;
                    int n2 = (inod+2)%4;
                    int n3 = (inod+3)%4;
                    prec* p0 = Vand[inod];
                    prec* p1 = Vand[n1]; prec* p2 = Vand[n2]; prec* p3 = Vand[n3];
                    //-------------------
                    temp_b = p1[1]*(p2[2]-p3[2]) + p2[1]*(p3[2]-p1[2]) + p3[1]*(p1[2] - p2[2]);
                    temp_b = -sign*temp_b;
                    //---
                    temp_c = p1[0]*(p2[2]-p3[2]) + p2[0]*(p3[2]-p1[2]) + p3[0]*(p1[2] - p2[2]);
                    temp_c = sign*temp_c;
                    //---
                    temp_d = p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1] - p2[1]);
                    temp_d = -sign*temp_d;
                    //--------------------
                    if (inod == 0) 
                    {
                        temp_a = p1[1]*(p2[2]*p3[0] - p2[0]*p3[2]) + p2[1]*(p3[2]*p1[0] - p3[0]*p1[2]) + p3[1]*(p1[2]*p2[0] - p1[0]*p2[2]);
                        temp_a = sign*temp_a;
                        detVand = temp_a + temp_b*p0[0] + temp_c*p0[1] + temp_d*p0[2];
                    }

                    Bloc[iel][inod] = temp_b/detVand;
                    Cloc[iel][inod] = temp_c/detVand;
                    Dloc[iel][inod] = temp_d/detVand;

                    sign = -sign;
                }
                Volume[iel] = std::abs(detVand)/6;
            }
            break;
        }
    }
    (*physics).Coef.resize(dim);
    switch (dim)
    {
    case 2:
        {
            (*physics).Coef[0].define((*physics).Bloc);
            (*physics).Coef[1].define((*physics).Cloc);
            break;
        }
    case 3:
        {
            (*physics).Coef[0].define((*physics).Bloc);
            (*physics).Coef[1].define((*physics).Cloc);
            (*physics).Coef[2].define((*physics).Dloc);
            break;
        }
    }
}



