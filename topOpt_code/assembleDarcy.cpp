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
    
    (*physics).Bloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Bloc = (*physics).Bloc.PP;
    (*physics).Cloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Cloc = (*physics).Cloc.PP;
    
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

void PROBLEM_DARCY::assemble() // i: usually refers to rows; j: usually refers to columns
{
    std::cout << "\n----------\n--| ASSEMBLE |--\n----------\n";
    // SYSMAT_DARCY::
    /*
       K + M/dT
    */

   // SYSMAT_DARCY BASE
   /*
       K
   */

    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    int nNodes = (*physics).nNodes;
    prec mu = (*physics).mu;
    std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    std::shared_ptr<int[]> elem_domain_id = (*physics).elem_geo_entities_ids.P;
    std::shared_ptr<prec*[]> Bloc = (*physics).Bloc.PP;
    std::shared_ptr<prec*[]> Cloc = (*physics).Cloc.PP;
    std::shared_ptr<prec*[]> Dloc = (*physics).Dloc.PP;

    nEvals = (dim+1)*(dim+1)*nElem;
    iSparse.initialize(nEvals);
    jSparse.initialize(nEvals);
    realPos.initialize(nEvals);
    
    std::shared_ptr<prec[]> coefH    = 0; 
    std::shared_ptr<prec[]> coefM    = 0;

    std::shared_ptr<int[]> i_sparse = VECTOR_INT::makePointer(nEvals);
    std::shared_ptr<int[]> j_sparse = VECTOR_INT::makePointer(nEvals);
    coefH     = VECTOR::makePointer(nEvals);
    coefM     = VECTOR::makePointer(nEvals);
    
    // switch(dim)
    // {
    //     case 2:  //  TWO DIMENSIONAL CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             2D                 
    //     {
    //         //---
    //         int count = 0;

    //         for (int iel = 0; iel < nElem; iel++) //enter the triangle
    //         {
    //             int temp_domain_id = elem_domain_id[iel];
    //             prec temp_permeability = domains_permeability[temp_domain_id];
    //             prec tempHFactor = Volume[iel] * temp_permeability / mu;

    //             int* tempElem = elem[iel];
    //             for (int iloc = 0; iloc < 3; iloc++) // select node i of the triangle
    //             {
    //                 int iglob = tempElem[iloc];
    //                 prec b_i = Bloc[iel][iloc];
    //                 prec c_i = Cloc[iel][iloc];
    //                 //----------------------------
    //                 for (int jloc = 0; jloc < 3; jloc++) // select node j of the triangle
    //                 {
    //                     int jglob = tempElem[jloc];
                        
    //                     // MASS term
    //                     if (iglob == jglob) coefM[count] = Volume[iel]/6;
    //                     else coefM[count] = Volume[iel]/12;

    //                     prec b_j = Bloc[iel][jloc];
    //                     prec c_j = Cloc[iel][jloc];

    //                     //  DIFFUSION term
    //                     coefH[count] = (b_i*b_j + c_i*c_j) * tempHFactor;
    //                     i_sparse[count] = iglob; 
    //                     j_sparse[count] = jglob; 
    //                     count++;
    //                 } // end jloc
    //             } // end iloc
    //         } // end iel
    //         //--------------------
    //         // BUILD SPARSE MATRIX
    //         //--------------------
    //         iSparse = i_sparse;
    //         jSparse = j_sparse;
            
    //         // build H matrix
    //         H.initialize(nNodes, nNodes, count, i_sparse, j_sparse, coefH);
    //         nTerms = H.nTerm;
    //         #pragma omp parallel for num_threads(PARALLEL::nThread)
    //         for (int ieval = 0; ieval < nEvals; ieval++)
    //         {
    //             realPos[ieval] = H.getPos(iSparse[ieval], jSparse[ieval]);
    //         }
    //         H.copyPattTo(M);

    //         // define M matrix
    //         M.defineZero(nNodes, nNodes);
    //         // #pragma omp parallel for num_threads(PARALLEL::nThread)
    //         for (int i = 0; i < nEvals; i++)
    //         {
    //             int iglob = i_sparse[i]; int jglob = j_sparse[i];
    //             M(iglob, jglob) += coefM[i];
    //         }
    //         break;
    //     }
    //     case 3:
    //     {
    //         //---
    //         int count = 0;
    //         for (int iel = 0; iel < nElem; iel++) //enter the triangle
    //         {
    //             int temp_domain_id = elem_domain_id[iel];
    //             prec temp_permeability = domains_permeability[temp_domain_id];
    //             prec tempHFactor = Volume[iel] * temp_permeability / mu;

    //             int* tempElem = elem[iel];
                
    //             for (int iloc = 0; iloc < 4; iloc++) // select node i of the tetra
    //             {
    //                 int iglob = tempElem[iloc];
    //                 prec b_i = Bloc[iel][iloc];
    //                 prec c_i = Cloc[iel][iloc];
    //                 prec d_i = Dloc[iel][iloc];
    //                 //----------------------------
    //                 for (int jloc = 0; jloc < 4; jloc++) // select node j of the tetra
    //                 {
    //                     int jglob = tempElem[jloc];
                        
    //                     // MASS term
    //                     if (iglob == jglob) coefM[count] = Volume[iel]/10;
    //                     else coefM[count] = Volume[iel]/20;

    //                     prec b_j = Bloc[iel][jloc];
    //                     prec c_j = Cloc[iel][jloc];
    //                     prec d_j = Dloc[iel][jloc];
    //                     // add DIFFUSION term
    //                     coefH[count] = (b_i*b_j + c_i*c_j + d_i*d_j) * tempHFactor;
    //                     i_sparse[count] = iglob; 
    //                     j_sparse[count] = jglob; 
    //                     count++;
    //                 } // end jloc
    //             } // end iloc
    //         } // end iel

    //         //--------------------
    //         // BUILD SPARSE MATRIX
    //         //--------------------
    //         iSparse = i_sparse;
    //         jSparse = j_sparse;
    //         H.initialize(nNodes, nNodes, count, i_sparse, j_sparse, coefH);
    //         nTerms = H.nTerm;
    //         #pragma omp parallel for num_threads(PARALLEL::nThread)
    //         for (int ieval = 0; ieval < nEvals; ieval++)
    //         {
    //             realPos[ieval] = H.getPos(iSparse[ieval], jSparse[ieval]);
    //         }
    //         H.copyPattTo(M);
    //         M.defineZero(nNodes, nNodes);
    //         // #pragma omp parallel for num_threads(PARALLEL::nThread)
    //         for (int i = 0; i < nEvals; i++)
    //         {
    //             int iglob = i_sparse[i]; int jglob = j_sparse[i];
    //             M(iglob, jglob) += coefM[i];
    //         }
    //         break;
    //     }
    // }
    
    //---
    int count = 0;
    for (int iel = 0; iel < nElem; iel++) //enter the triangle
    {
        int* tempElem = elem[iel];
        int temp_domain_id = elem_domain_id[iel];
        prec temp_permeability = domains_permeability[temp_domain_id];
        prec tempHFactor = Volume[iel] / mu * temp_permeability;
        prec tempMFactor = Volume[iel]/ ((dim+1)*(dim+2));

        //----------------------------
        for (int iloc = 0; iloc < dim+1; iloc++) // select node i of the tetra
        {
            int iglob = tempElem[iloc];
            VECTOR i_coef(dim);
            for (int icoef = 0; icoef < dim; icoef++)
            {
                i_coef[icoef] = (*physics).Coef[icoef][iel][iloc];
            }

            //----------------------------
            for (int jloc = 0; jloc < dim+1; jloc++) // select node j of the tetra
            {
                int jglob = tempElem[jloc];

                // MASS term
                if (iglob == jglob) coefM[count] = 2 * tempMFactor;
                else coefM[count] = tempMFactor;
                
                // DIFFUSION term
                VECTOR j_coef(dim);
                for (int jcoef = 0; jcoef < dim; jcoef++)
                {
                    j_coef[jcoef] = (*physics).Coef[jcoef][iel][jloc];
                }
                prec tempHCoef = i_coef.dot(j_coef);
                coefH[count] = tempHCoef * tempHFactor;
                i_sparse[count] = iglob; 
                j_sparse[count] = jglob; 
                count++;
            } // end jloc
        } // end iloc
    } // end iel

    //--------------------
    // BUILD SPARSE MATRIX
    //--------------------
    iSparse = i_sparse;
    jSparse = j_sparse;
    H.initialize(nNodes, nNodes, count, i_sparse, j_sparse, coefH);
    nTerms = H.nTerm;
    #pragma omp parallel for num_threads(PARALLEL::nThread)
    for (int ieval = 0; ieval < nEvals; ieval++)
    {
        realPos[ieval] = H.getPos(iSparse[ieval], jSparse[ieval]);
    }
    H.copyPattTo(M);
    M.defineZero(nNodes, nNodes);
    for (int i = 0; i < nEvals; i++)
    {
        int iglob = i_sparse[i]; int jglob = j_sparse[i];
        M(iglob, jglob) += coefM[i];
    }

    SYSMAT_base = H;
}

void PROBLEM_DARCY::addToSysmat(CSRMAT &mat, prec factor)
{
    bool same_path = SYSMAT.checkPattern(mat);
    if (same_path)
    {
        std::shared_ptr<prec[]> coefM = M.coef;
        for (int iterm = 0; iterm < SYSMAT.nTerm; iterm++)
        {
            SYSMAT.coef[iterm] += coefM[iterm] * factor;
        }
    }
    else
    {
        std::cout << "\n\n !!! SYSMAT and MASS MATRIX have different patterns !!! \n\n";
        throw_line("ERROR: not handled SYSMAT sum case\n");
    }
}


