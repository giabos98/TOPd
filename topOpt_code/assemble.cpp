#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"

void PROBLEM_NS::localBasis()
{
    std::cout << "\n----------\n--| LOCAL BASIS  |--\n----------\n";
    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    std::shared_ptr<prec*[]>  coord = (*physics).coord.PP; 
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    (*physics).Volume.initialize(nElem);
    std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    (*physics).Volume_v.initialize(nElem_v);
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    
    // PRESSURE
    (*physics).Bloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Bloc = (*physics).Bloc.PP;
    (*physics).Cloc.initialize(nElem, dim+1);
    std::shared_ptr<prec*[]> Cloc = (*physics).Cloc.PP;
    // VELOCITY
    (*physics).Bloc_v.initialize(nElem_v, dim+1);
    std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    (*physics).Cloc_v.initialize(nElem_v, dim+1);
    std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    

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
            // VELOCITY NODES
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* locNodes = elem_v[iel];
                int node0 = locNodes[0]; int node1 = locNodes[1]; int node2 = locNodes[2];
                prec p0x = coord_v[node0][0]; prec p0y = coord_v[node0][1];
                prec p1x = coord_v[node1][0]; prec p1y = coord_v[node1][1];
                prec p2x = coord_v[node2][0]; prec p2y = coord_v[node2][1];

                prec detVand = p2y*(p1x-p0x) - p1y*(p2x-p0x) + p0y*(p2x-p1x);
                Volume_v[iel] = std::abs(detVand)/2;
                for (int inod = 0; inod < 3; inod++)
                {
                    int n1 = (inod+1)%3;
                    int n2 = (inod+2)%3;
                    Bloc_v[iel][inod] = (coord_v[locNodes[n1]][1] - coord_v[locNodes[n2]][1])/detVand;
                    Cloc_v[iel][inod] = (coord_v[locNodes[n2]][0] - coord_v[locNodes[n1]][0])/detVand;
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
            // VELOCITY
            (*physics).Dloc_v.initialize(nElem_v, 4);
            std::shared_ptr<prec*[]> Dloc_v = (*physics).Dloc_v.PP;
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* locNodes = elem_v[iel];
                MATRIX Vand(4,dim);
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        Vand[i][j] = coord_v[locNodes[i]][j];
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

                    Bloc_v[iel][inod] = temp_b/detVand;
                    Cloc_v[iel][inod] = temp_c/detVand;
                    Dloc_v[iel][inod] = temp_d/detVand;

                    sign = -sign;
                }
                Volume_v[iel] = std::abs(detVand)/6;
            }
            break;
        }
    }
    (*physics).Coef_v.resize(dim);
    switch (dim)
    {
    case 2:
        {
            (*physics).Coef_v[0].define((*physics).Bloc_v);
            (*physics).Coef_v[1].define((*physics).Cloc_v);
            break;
        }
    case 3:
        {
            (*physics).Coef_v[0].define((*physics).Bloc_v);
            (*physics).Coef_v[1].define((*physics).Cloc_v);
            (*physics).Coef_v[2].define((*physics).Dloc_v);
            break;
        }
    }
}


void PROBLEM_NS::assemble() // i: usually refers to rows; j: usually refers to columns
{
    std::cout << "\n----------\n--| ASSEMBLE |--\n----------\n";
    // SYSMAT_NS::
    /*
       K + N(u) + M/dT,                0,                     0,      B1^T

                     0,  K + N(u) + M/dT,                     0,      B2^T

                     0,                 0,      K + N(u) + M/dT,      B3^T

                    B1,                B2,                   B3,         0
    */

   // SYSMAT_NS BASE
   /*
       K,     0,    0,   B1^T

       0,     K,    0,   B2^T

       0,     0,    K,   B3^T

      B1,    B2,   B3,      0
   */

    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    int nElem_v = (*physics).nElem_v;
    int nNodes = (*physics).nNodes;
    int nNodes_v = (*physics).nNodes_v;
    prec rho = (*physics).rho;
    prec ni = (*physics).ni;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    std::shared_ptr<prec*[]> Dloc_v = (*physics).Dloc_v.PP;

    nEvals = (dim+1)*(dim+1)*nElem_v;
    iSparse.initialize(nEvals);
    jSparse.initialize(nEvals);
    realPos.initialize(nEvals);
    int nEvalsB = 4*(dim-1)*(dim+1)*(dim+1)*nElem;
    
    std::shared_ptr<prec[]> coefH    = 0; 
    std::shared_ptr<prec[]> coefM    = 0;

    std::shared_ptr<int[]> i_sparse = VECTOR_INT::makePointer(nEvals);
    std::shared_ptr<int[]> j_sparse = VECTOR_INT::makePointer(nEvals);
    coefH     = VECTOR::makePointer(nEvals);
    coefM     = VECTOR::makePointer(nEvals);
    
    switch(dim)
    {
    case 2:  //  TWO DIMENSIONAL CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             2D                 
    {
        //---
        int count = 0;

        for (int iel = 0; iel < nElem_v; iel++) //enter the triangle
        {
            prec tempHFactor = Volume_v[iel] * ni;

            int* tempElem = elem_v[iel];
            for (int iloc = 0; iloc < 3; iloc++) // select node i of the triangle
            {
                int iglob = tempElem[iloc];
                prec b_i = Bloc_v[iel][iloc];
                prec c_i = Cloc_v[iel][iloc];
                //----------------------------
                for (int jloc = 0; jloc < 3; jloc++) // select node j of the triangle
                {
                    int jglob = tempElem[jloc];
                    
                    // MASS term
                    if (iglob == jglob) coefM[count] = Volume_v[iel]/6;
                    else coefM[count] = Volume_v[iel]/12;

                    prec b_j = Bloc_v[iel][jloc];
                    prec c_j = Cloc_v[iel][jloc];
                    // add DIFFUSION term
                    coefH[count] = (b_i*b_j + c_i*c_j) * tempHFactor;
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
        
        // build H matrix
        H.initialize(nNodes_v, nNodes_v, count, i_sparse, j_sparse, coefH);
        nTerms = H.nTerm;
        #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int ieval = 0; ieval < nEvals; ieval++)
        {
            realPos[ieval] = H.getPos(iSparse[ieval], jSparse[ieval]);
        }
        H.copyPattTo(M);

        // define M matrix
        M.defineZero(nNodes_v, nNodes_v);
        // #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int i = 0; i < nEvals; i++)
        {
            int iglob = i_sparse[i]; int jglob = j_sparse[i];
            M(iglob, jglob) += coefM[i];
        }

        // initialize N matrix coef
        //--------------------------------------------------------------------------
        // BUILD PRESSURE GRADIENT MATRIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //--------------------------------------------------------------------------
        i_sparse.reset();
        i_sparse = VECTOR_INT::makePointer(nEvalsB);
        std::shared_ptr<int[]> k_sparse = VECTOR_INT::makePointer(nEvalsB);
        std::shared_ptr<prec[]> coefB1 = VECTOR::makePointer(nEvalsB);
        std::shared_ptr<prec[]> coefB2 = VECTOR::makePointer(nEvalsB);
        
        count = 0;
        for (int iel = 0; iel < nElem; iel++)
        {
            int* tempPressureElem = elem[iel];
            for (int kloc = 0; kloc < 3; kloc++)
            {
                int kglob = tempPressureElem[kloc];
                for (int iel_v = 0; iel_v < 4; iel_v++)
                {
                    int globElem_v = 4*iel + iel_v;
                    int* tempVelElem = elem_v[globElem_v];
                    prec tempFactor = Volume_v[globElem_v]/3;
                    for (int iloc = 0; iloc < 3; iloc++)
                    {
                        int iglob = tempVelElem[iloc];
                        prec b_i = Bloc_v[globElem_v][iloc];
                        prec c_i = Cloc_v[globElem_v][iloc];
                        prec tempValueB1 = 0;
                        prec tempValueB2 = 0;
                        //-----------------------------
                        if (kloc == iel_v)
                        {
                            tempValueB1 += b_i * tempFactor;
                            tempValueB2 += c_i * tempFactor;
                        }
                        //-----------------------------
                        if (kloc == iel_v || iel_v == 3)
                        {
                            tempValueB1 += b_i * tempFactor;
                            tempValueB2 += c_i * tempFactor;
                        } else
                        {
                            tempValueB1 += 0.5 * b_i * tempFactor;
                            tempValueB2 += 0.5 * c_i * tempFactor;
                        }
                        // update
                        i_sparse[count] = iglob;
                        k_sparse[count] = kglob;
                        coefB1[count] = -tempValueB1/rho;
                        coefB2[count] = -tempValueB2/rho;
                        count++;
                    }
                }
            }
        }
        //------------------------------------------------------------
        CSRMAT B1(nNodes, nNodes_v, count, k_sparse, i_sparse, coefB1);
        CSRMAT B2;

        nTermB = B1.nTerm;;
        B1.copyPattTo(B2);
        B2.defineZero(nNodes, nNodes_v);
        // #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int i = 0; i < nEvals; i++)
        {
            int iglob = i_sparse[i]; int kglob = k_sparse[i];
            B2(kglob, iglob) += coefB2[i];
        }

        //-------------------------------------------------------------
        // ASSEMBLY OF SYSMAT_NS BASE
        //-------------------------
        B.resize(dim);
        B[0] = B1; B[1] = B2;

        createSysmatBase(H, B);
        break;
    }
    case 3:
    {
        //---
        int count = 0;
        for (int iel = 0; iel < nElem_v; iel++) //enter the triangle
        {
            prec tempHFactor = Volume_v[iel] * ni;

            int* tempElem = elem_v[iel];
            
            for (int iloc = 0; iloc < 4; iloc++) // select node i of the tetra
            {
                int iglob = tempElem[iloc];
                prec b_i = Bloc_v[iel][iloc];
                prec c_i = Cloc_v[iel][iloc];
                prec d_i = Dloc_v[iel][iloc];
                //----------------------------
                for (int jloc = 0; jloc < 4; jloc++) // select node j of the tetra
                {
                    int jglob = tempElem[jloc];
                    
                    // MASS term
                    if (iglob == jglob) coefM[count] = Volume_v[iel]/10;
                    else coefM[count] = Volume_v[iel]/20;

                    prec b_j = Bloc_v[iel][jloc];
                    prec c_j = Cloc_v[iel][jloc];
                    prec d_j = Dloc_v[iel][jloc];
                    // add DIFFUSION term
                    coefH[count] = (b_i*b_j + c_i*c_j + d_i*d_j) * tempHFactor;
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
        H.initialize(nNodes_v, nNodes_v, count, i_sparse, j_sparse, coefH);
        nTerms = H.nTerm;
        #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int ieval = 0; ieval < nEvals; ieval++)
        {
            realPos[ieval] = H.getPos(iSparse[ieval], jSparse[ieval]);
        }
        H.copyPattTo(M);
        M.defineZero(nNodes_v, nNodes_v);
        // #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int i = 0; i < nEvals; i++)
        {
            int iglob = i_sparse[i]; int jglob = j_sparse[i];
            M(iglob, jglob) += coefM[i];
        }


        //--------------------------------------------------------------------------
        // BUILD PRESSURE GRADIENT MATRIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //--------------------------------------------------------------------------
        i_sparse.reset();
        i_sparse = VECTOR_INT::makePointer(nEvalsB);
        std::shared_ptr<int[]> k_sparse = VECTOR_INT::makePointer(nEvalsB);
        std::shared_ptr<prec[]> coefB1 = VECTOR::makePointer(nEvalsB);
        std::shared_ptr<prec[]> coefB2 = VECTOR::makePointer(nEvalsB);
        std::shared_ptr<prec[]> coefB3 = VECTOR::makePointer(nEvalsB);

        VECTOR factors(4); factors[0] = 3.0; factors[1] = 2.0; factors[2] = 1.0; factors[3] = 2.0;
        count = 0;
        for (int iel = 0; iel < nElem; iel++)
        {
            int* tempPressureElem = elem[iel];
            for (int kloc = 0; kloc < 4; kloc++)
            {
                int kglob = tempPressureElem[kloc];
                for (int iel_v = 0; iel_v < 8; iel_v++)
                {
                    int globElem_v = 8*iel + iel_v;
                    int* tempVelElem = elem_v[globElem_v];
                    prec tempFactor = Volume_v[globElem_v]/4;

                    for (int iloc = 0; iloc < 4; iloc++)
                    {
                        int iglob = tempVelElem[iloc];
                        prec b_i = Bloc_v[globElem_v][iloc];
                        prec c_i = Cloc_v[globElem_v][iloc];
                        prec d_i = Dloc_v[globElem_v][iloc];
                        prec tempValueB1;
                        prec tempValueB2;
                        prec tempValueB3;
                        //-----------------------------
                        if (kloc == iel_v)
                        {
                            tempValueB1 = 2.5 * b_i * tempFactor; // (1.0 + 0.5 * 3) = 2.5
                            tempValueB2 = 2.5 * c_i * tempFactor; // (1.0 + 0.5 * 3) = 2.5
                            tempValueB3 = 2.5 * d_i * tempFactor; // (1.0 + 0.5 * 3) = 2.5
                        }
                        else if (iel_v < 4)
                        {
                            tempValueB1 = 0.5 * b_i * tempFactor;
                            tempValueB2 = 0.5 * c_i * tempFactor;
                            tempValueB3 = 0.5 * d_i * tempFactor;
                        }
                        else
                        {
                            int idFactor = (iel_v - kloc) % 4;
                            prec factor = factors[idFactor];
                            tempValueB1 = 0.5 * factor * b_i * tempFactor;
                            tempValueB2 = 0.5 * factor * c_i * tempFactor;
                            tempValueB3 = 0.5 * factor * d_i * tempFactor;
                        }
                        // update
                        i_sparse[count] = iglob;
                        k_sparse[count] = kglob;
                        coefB1[count] = -tempValueB1/rho;
                        coefB2[count] = -tempValueB2/rho;
                        coefB3[count] = -tempValueB3/rho;
                        count++;
                    }
                }
            }
        }
        //------------------------------------------------------------
        CSRMAT B1(nNodes, nNodes_v, count, k_sparse, i_sparse, coefB1);
        CSRMAT B2; CSRMAT B3;

        nTermB = B1.nTerm;
        B1.copyPattTo(B2); B1.copyPattTo(B3);
        B2.defineZero(nNodes, nNodes_v); B3.defineZero(nNodes, nNodes_v);
        // #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int i = 0; i < nEvalsB; i++)
        {
            int iglob = i_sparse[i]; int kglob = k_sparse[i];
            B2(kglob, iglob) += coefB2[i];
            B3(kglob, iglob) += coefB3[i];
        }

        //-------------------------------------------------------------
        // ASSEMBLY OF SYSMAT_NS BASE
        //-------------------------
        B.resize(dim);
        B[0] = B1; B[1] = B2; B[2] = B3;


        createSysmatBase(H, B);
        B1.dlt(); B2.dlt(); B3.dlt();
        break;
    }
    }

    //initialize other matrices
    N.initialize(nTerms);
    Ma.initialize(nTerms);
}

//-----------------------------------------------------------------------------
// SYSMAT_NS BASE
//-----------------------------------------------------------------------------
void PROBLEM_NS::createSysmatBase(CSRMAT &H, std::vector<CSRMAT> &B)
{
    int dim = (*physics).dim;
    int nNodes = (*physics).nNodes;
    int nNodes_v = (*physics).nNodes_v;

    int nTerm = dim * (H.nTerm + 2*B[0].nTerm);
    int BLastRow = nNodes;
    if (nNeuBound == 0)
    {
        BLastRow = nNodes-1;
    }

    realPosSysmat.resize(1);
    realPosSysmat[0].resize(dim);

    for (int idim = 0; idim < dim; idim++) realPosSysmat[0][idim].initialize(nTerms);

    int n = nNodes_v*dim + nNodes;
    SYSMAT_base.initialize(n, n, nTerm);
    std::shared_ptr<int[]> iat    = SYSMAT_base.iat;
    iat[0] = 0;
    std::shared_ptr<int[]> ja     = SYSMAT_base.ja;
    std::shared_ptr<prec[]> coef  = SYSMAT_base.coef;

    std::shared_ptr<int[]> iatH  = H.iat;
    std::shared_ptr<int[]> jaH   = H.ja;
    std::shared_ptr<prec[]> coefH = H.coef;

    int pos = 0;
    int globRow = 0;
    for (int icomp = 0; icomp < dim ; icomp++)
    {
        CSRMAT tempB = B[icomp].getTranspose();
        std::shared_ptr<int[]> iatB = tempB.iat;
        std::shared_ptr<int[]> jaB = tempB.ja;
        std::shared_ptr<prec[]> coefB = tempB.coef;
        for (int irow = 0; irow < nNodes_v; irow++)
        {
            for (int ipos_v = iatH[irow]; ipos_v < iatH[irow+1]; ipos_v++)
            {
                prec tempPos = jaH[ipos_v] + icomp*nNodes_v;
                ja[pos] = tempPos;
                coef[pos] = coefH[ipos_v];
                pos++;
            }
            //----
            for (int ipos_p = iatB[irow]; ipos_p < iatB[irow+1]; ipos_p++)
            {
                ja[pos] = jaB[ipos_p] + dim*nNodes_v;
                coef[pos] = coefB[ipos_p];
                pos++;
            }
            globRow++;
            iat[globRow] = pos;
        }
    }
    //-----------------------------------------------------------------------
    for (int irow = 0; irow < BLastRow; irow++)
    {
        for (int icomp = 0; icomp < dim; icomp++)
        {
            CSRMAT tempB = B[icomp];
            std::shared_ptr<int[]> iatB = tempB.iat;
            std::shared_ptr<int[]> jaB = tempB.ja;
            std::shared_ptr<prec[]> coefB = tempB.coef;
            for (int ipos = iatB[irow]; ipos < iatB[irow+1]; ipos++)
            {
                ja[pos] = jaB[ipos] + icomp*nNodes_v;
                coef[pos] = coefB[ipos];
                pos++;
            }
        }
        globRow++;
        iat[globRow] = pos;
    }
    //----------------------
    // INSERT ZERO MEAN PRESSURE CONDITION IF NO NEU BOUND
    //----------------------
    if (nNeuBound == 0)
    {
        for (int icomp = 0; icomp < dim; icomp++)
        {
            CSRMAT tempB = B[icomp];
            std::shared_ptr<int[]> iatB = tempB.iat;
            std::shared_ptr<int[]> jaB = tempB.ja;
            std::shared_ptr<prec[]> coefB = tempB.coef;
            for (int ipos = iatB[BLastRow]; ipos < iatB[BLastRow+1]; ipos++)
            {
                ja[pos] = jaB[ipos] + icomp*nNodes_v;
                coef[pos] = coefB[ipos];;
                pos++;
            }
        }
        globRow++;
        iat[globRow] = pos;
    }

    if (nSymm > 0) 
    {
        CSRMAT tempMat = SYSMAT_base + J;
        SYSMAT_base = tempMat;
    }

    // COMPUTE realPosSYSMAT
    iat = SYSMAT_base.iat; ja = SYSMAT_base.ja;
    for (int icomp = 0; icomp < dim; icomp++)
    {
        int ipos_v = 0;
        for (int irow = icomp*nNodes_v; irow < (icomp+1)*nNodes_v; irow++)
        {
            int startPos = iat[irow]; int endPos = iat[irow+1];

            for (int pos = startPos; pos < endPos; pos++)
            {
                int tempCol = ja[pos];
                if (tempCol < nNodes_v*icomp) continue;
                else if (tempCol >= nNodes_v*icomp + nNodes_v) break;
                realPosSysmat[0][icomp][ipos_v] = pos;
                ipos_v++;
            }
        }
    }
                    
}

//--------------------------------------------------------
// BUILD CONVECTIVE MATRIX N
//--------------------------------------------------------
void PROBLEM_NS::assembleN(VECTOR &lastSol, CSRMAT& SYSMAT_final, VECTOR &rhs_final) //(CSRMAT &SYSMAT_final, VECTOR &rhs_final)
{
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    prec ni = (*physics).ni;
    std::shared_ptr<prec[]> h_v = (*physics).h_v.P;
    std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    std::shared_ptr<prec*[]> Dloc_v = (*physics).Dloc_v.PP;


    //----------------------------------------------------------------------
    // BUILD N+S
    //----------------------------------------------------------------------

    std::shared_ptr<prec*[]> U_r = VECTOR::makeDoublePointer(dim);

    for (int icomp = 0; icomp < dim; icomp++)
    {
        U_r[icomp] = &lastSol[icomp*nNodes_v];
    }
    prec* u0_r = U_r[0];
    prec* u1_r = U_r[1];
    prec* u2_r = nullptr;

    //prec* coef = SYSMAT_final.coef;

    //int nEvals = (dim+1) * (dim+1) * nElem_v;
 
    std::shared_ptr<int[]> i_sparse  = VECTOR_INT::makePointer(nEvals);
    std::shared_ptr<int[]> j_sparse  = VECTOR_INT::makePointer(nEvals);
    std::shared_ptr<prec[]> coefN    = VECTOR::makePointer(nEvals);

    

    VECTOR tau(nElem_v);
    int count = 0;
    switch (dim)
    {
        case 2:
        {
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* tempElem = elem_v[iel];
                prec tempFactor = Volume_v[iel] / 3;

                VECTOR beta0(3);
                VECTOR beta1(3);
                VECTOR beta0_mean(3); 
                VECTOR beta1_mean(3); 
                prec sum0 = 0;
                prec sum1 = 0;
                for (int iloc = 0; iloc < 3; iloc++)
                {
                    int iglob = tempElem[iloc];
                    prec tempBeta0 = u0_r[iglob];
                    prec tempBeta1 = u1_r[iglob];
                    
                    sum0 += tempBeta0; 
                    sum1 += tempBeta1;

                    beta0[iloc] = tempBeta0;
                    beta1[iloc] = tempBeta1;
                    beta0_mean[iloc] = tempBeta0;
                    beta1_mean[iloc] = tempBeta1;
                }
                beta0_mean += sum0;
                beta1_mean += sum1;

                beta0_mean /= 4;
                beta1_mean /= 4;


                prec tempBetaNorm = sqrt(sum0*sum0 + sum1*sum1)/3;

                prec tempTau;
                if (tempBetaNorm < 1e-8) tempTau = 0;
                else 
                {
                    /*
                    Pe_h = beta*h/ni_tot ---> ni_tot > beta*h/1.8

                    ni_tot = beta^2*tau + ni ---> tau = (h*beta/1.8-ni)/beta/beta

                    
                    */
                    // tempTau = 0.01/h_v[iel] /ni/2/tempBetaNorm;
                    tempTau = (h_v[iel]/(1.8*tempBetaNorm) - ni/(tempBetaNorm*tempBetaNorm));    
                    if (tempTau < 0) tempTau = 0;
                    else if (tempTau > deltaT/2) tempTau = deltaT/2;
                }
                if (tempTau > 0.01) 
                {
                    tempTau = 0.01;
                }
                tau[iel] = tempTau; 
                //------------------------------

                for (int iloc = 0; iloc < 3; iloc++)
                {
                    //------------------------------
                    // SELECT FIRST VERTEX
                    //------------------------------
                    int iglob = tempElem[iloc];

                    prec b_i = Bloc_v[iel][iloc];
                    prec c_i = Cloc_v[iel][iloc];
                    for (int jloc = 0; jloc < 3; jloc++)
                    {
                        int jglob = tempElem[jloc];
                        prec b_j = Bloc_v[iel][jloc];
                        prec c_j = Cloc_v[iel][jloc];
                        //--------------------------
                        prec tempBeta0_mean = beta0_mean[iloc];
                        prec tempBeta1_mean = beta1_mean[iloc];
                        prec Nval = (tempBeta0_mean * b_j + tempBeta1_mean * c_j) * tempFactor;

                        // ADD SUPG STABILIZATION
                        prec Sval = 0;

                        //  mass term
                        // Sval += (b_i*beta0_mean[jloc] + c_i*beta1_mean[jloc])/ deltaT;   

                        // numerical diffusion term
                        Sval += (b_j*b_i*beta0.dot(beta0_mean) + (b_j*c_i + c_j*b_i)*beta0.dot(beta1_mean) + c_i*c_j*beta1.dot(beta1_mean)); 

                        Sval *= tempFactor; Sval *= tempTau;
                        i_sparse[count] = iglob;
                        j_sparse[count] = jglob;
                        coefN[count]     = Nval + Sval;
                        count++;
                    }
                }
            }
            break;
        }
        case 3:
        {
            u2_r = U_r[2];
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* tempElem = elem_v[iel];
                prec tempFactor = Volume_v[iel] / 4;

                VECTOR beta0(4);
                VECTOR beta1(4);
                VECTOR beta2(4);
                VECTOR beta0_mean(4); 
                VECTOR beta1_mean(4); 
                VECTOR beta2_mean(4); 
                prec sum0 = 0;
                prec sum1 = 0;
                prec sum2 = 0;
                for (int iloc = 0; iloc < 4; iloc++)
                {
                    int iglob = tempElem[iloc];
                    prec tempBeta0 = u0_r[iglob];
                    prec tempBeta1 = u1_r[iglob];
                    prec tempBeta2 = u2_r[iglob];
                    
                    sum0 += tempBeta0; 
                    sum1 += tempBeta1;
                    sum2 += tempBeta2;

                    beta0[iloc] = tempBeta0;
                    beta1[iloc] = tempBeta1;
                    beta2[iloc] = tempBeta2;
                    beta0_mean[iloc] = tempBeta0;
                    beta1_mean[iloc] = tempBeta1;
                    beta2_mean[iloc] = tempBeta2;
                }
                beta0_mean += sum0;
                beta1_mean += sum1;
                beta2_mean += sum2;

                beta0_mean /= 5.0;
                beta1_mean /= 5.0;
                beta2_mean /= 5.0;

                


                prec tempBetaNorm = sqrt(sum0*sum0 + sum1*sum1 + sum2*sum2)/4; //NORM OF THE MEAN

                prec tempTau;
                if (tempBetaNorm < 1e-8) tempTau = 0;
                else 
                {
                    /*
                    Pe_h = beta*h/ni_tot ---> ni_tot > beta*h/1.8

                    ni_tot = beta^2*tau + ni ---> tau = (h*beta/1.8-ni)/beta/beta

                    
                    */
                    // tempTau = 0.01/h_v[iel] /ni/2/tempBetaNorm;
                    tempTau = (h_v[iel]/(1.8*tempBetaNorm) - ni/(tempBetaNorm*tempBetaNorm));    
                    if (tempTau < 0) tempTau = 0;
                    else if (tempTau > deltaT/2) tempTau = deltaT/2;
                }
                if (tempTau > 0.01) 
                {
                    tempTau = 0.01;
                }
                tau[iel] = tempTau; 
                //------------------------------

                for (int iloc = 0; iloc < 4; iloc++)
                {
                    //------------------------------
                    // SELECT FIRST VERTEX
                    //------------------------------
                    int iglob = tempElem[iloc];

                    prec b_i = Bloc_v[iel][iloc];
                    prec c_i = Cloc_v[iel][iloc];
                    prec d_i = Dloc_v[iel][iloc];
                    for (int jloc = 0; jloc < 4; jloc++)
                    {
                        int jglob = tempElem[jloc];
                        prec b_j = Bloc_v[iel][jloc];
                        prec c_j = Cloc_v[iel][jloc];
                        prec d_j = Dloc_v[iel][jloc];
                        //--------------------------
                        prec tempBeta0_mean = beta0_mean[iloc];
                        prec tempBeta1_mean = beta1_mean[iloc];
                        prec tempBeta2_mean = beta2_mean[iloc];
                        prec Nval = (tempBeta0_mean * b_j + tempBeta1_mean * c_j + tempBeta2_mean * d_j) * tempFactor;

                        // ADD SUPG STABILIZATION
                        prec Sval = 0;

                        //  mass term
                        // Sval += (b_i*beta0_mean[jloc] + c_i*beta1_mean[jloc] + d_i*beta2_mean[jloc])/ deltaT;   

                        // numerical diffusion term
                        Sval += b_j*b_i*beta0.dot(beta0_mean) + c_i*c_j*beta1.dot(beta1_mean) + d_i*d_j*beta2.dot(beta2_mean); 

                        Sval += (b_j*c_i + c_j*b_i)*beta0.dot(beta1_mean) + (c_j*d_i + d_j*c_i)*beta1.dot(beta2_mean) + (d_j*b_i + b_j*d_i)*beta2.dot(beta0_mean);

                        Sval *= tempFactor; Sval *= tempTau;
                        i_sparse[count] = iglob;
                        j_sparse[count] = jglob;
                        coefN[count]    = Nval + Sval;
                        count++;
                    }
                }
            }
            break;
        }
    }
    // for (int locPos = 0; locPos < count; locPos++)
    // {
    //     int iglob = i_sparse[locPos];
    //     int jglob = j_sparse[locPos];
    //     prec Nval = N[locPos];
    //     int posGlob = SYSMAT_final.getPos(iglob, jglob);

    //     for (int icomp = 0; icomp < dim; icomp++)
    //     {
    //         int tempPosGlob = posGlob + rescaleNS*icomp;
    //         coef[tempPosGlob] += Nval;
    //     }
    // }
    // std::cout << nTerms << "\n";
    // N.initialize(nTerms);
    // std::cout << N.length << "\n";

    resizeCoef(coefN, N);
 

    if (tau.norm() < 1e-9)
    {
        return;
    }
    //--------------------------------------------------------------
    // build S_B
    //--------------------------------------------------------------

    // std::shared_ptr<prec[]> coef = SYSMAT_final.coef;
    // int rescaleCoef = (SYSMAT_final.nTerm-dim*H.nTerm)/2;
    // int countB = 0;
    // int nElem = (*physics).nElem;
    // MATRIX_INT elem = (*physics).elem;
    // MATRIX Bloc = (*physics).Bloc;
    // MATRIX Cloc = (*physics).Cloc;
    // MATRIX Dloc = (*physics).Dloc;
    // prec rho = (*physics).rho;
    // switch (dim)
    // {
    //     case 2:
    //     {
    //         int nEvals = nElem*3*4*3;
    //         std::shared_ptr<int []> i_sparse(new int[nEvals]);
    //         std::shared_ptr<int []> k_sparse(new int[nEvals]);
    //         std::shared_ptr<prec []> coefSB1(new prec[nEvals]);
    //         std::shared_ptr<prec []> coefSB2(new prec[nEvals]);
            
    //         for (int iel = 0; iel < nElem; iel++)
    //         {
    //             int* tempElem = elem[iel];

    //             prec* tempBloc = Bloc[iel];
    //             prec* tempCloc = Cloc[iel];
    //             for (int kloc = 0; kloc < 3; kloc++)
    //             {
    //                 int kglob = tempElem[kloc];
    //                 prec b_k = tempBloc[kloc];
    //                 prec c_k = tempCloc[kloc];
    //                 for (int iel_v = 0; iel_v < 4; iel_v++)
    //                 {

    //                     int globElemId = iel*4+iel_v;
    //                     int* tempElem_v = elem_v[globElemId];
    //                     //---
    //                     prec* tempBloc_v = Bloc_v[globElemId];
    //                     prec* tempCloc_v = Cloc_v[globElemId];

    //                     prec tempFactor = Volume_v[globElemId]/(3*rho) * tau[globElemId];
    //                     //---
    //                     for (int iloc = 0; iloc < 3; iloc++)
    //                     {
    //                         prec b_i = tempBloc_v[iloc];
    //                         prec c_i = tempCloc_v[iloc];
    //                         //---
    //                         int iglob = tempElem_v[iloc];
    //                         prec tempU_0 = u0_r[iglob];
    //                         prec tempU_1 = u1_r[iglob];

    //                         coefSB1[countB] = b_k*(b_i*tempU_0 + c_i*tempU_1)*tempFactor;
    //                         coefSB2[countB] = c_k*(b_i*tempU_0 + c_i*tempU_1)*tempFactor;
    //                         i_sparse[countB] = iglob;
    //                         k_sparse[countB] = kglob;
    //                         countB++;
    //                     }
    //                 }
    //             }
    //         }
    //         //----------------------------------
    //         // ADD SB1 and SB2 to SYSMAT_final          

    //         // VECTOR::print2(coefSB1, coefSB2, countB, "coefSB1", "coefSB2");

    //         for (int locPos = 0; locPos < countB; locPos++)
    //         {
    //             int iglob = i_sparse[locPos];
    //             int kglob = k_sparse[locPos] + dim*nNodes_v;
    //             int pos = SYSMAT_final.getPos(iglob, kglob);
    //             coef[pos] += coefSB1[pos];
    //             coef[pos+rescaleCoef] += coefSB2[pos];
    //         }
    //         break;
    //     }
    //     case 3:
    //     {
    //         int nEvals = nElem*4*8*4;
    //         std::shared_ptr<int []> i_sparse(new int[nEvals]);
    //         std::shared_ptr<int []> k_sparse(new int[nEvals]);
    //         std::shared_ptr<prec []> coefSB1(new prec[nEvals]);
    //         std::shared_ptr<prec []> coefSB2(new prec[nEvals]);
    //         std::shared_ptr<prec []> coefSB3(new prec[nEvals]);
    //         for (int iel = 0; iel < nElem; iel++)
    //         {
    //             int* tempElem = elem[iel];

    //             prec* tempBloc = Bloc[iel];
    //             prec* tempCloc = Cloc[iel];
    //             prec* tempDloc = Dloc[iel];
    //             for (int kloc = 0; kloc < 4; kloc++)
    //             {
    //                 int kglob = tempElem[kloc];
    //                 prec b_k = tempBloc[kloc];
    //                 prec c_k = tempCloc[kloc];
    //                 prec d_k = tempDloc[kloc];
    //                 for (int iel_v = 0; iel_v < 8; iel_v++)
    //                 {

    //                     int globElemId = iel*8+iel_v;
    //                     int* tempElem_v = elem_v[globElemId];
    //                     //---
    //                     prec* tempBloc_v = Bloc_v[globElemId];
    //                     prec* tempCloc_v = Cloc_v[globElemId];
    //                     prec* tempDloc_v = Dloc_v[globElemId];

    //                     prec tempFactor = Volume_v[globElemId]/(4*rho) * tau[globElemId];
    //                     //---
    //                     for (int iloc = 0; iloc < 4; iloc++)
    //                     {
    //                         prec b_i = tempBloc_v[iloc];
    //                         prec c_i = tempCloc_v[iloc];
    //                         prec d_i = tempDloc_v[iloc];
    //                         //---
    //                         int iglob = tempElem_v[iloc];
    //                         prec tempU_0 = u0_r[iglob];
    //                         prec tempU_1 = u1_r[iglob];
    //                         prec tempU_2 = u2_r[iglob];

    //                         coefSB1[countB] = b_k*(b_i*tempU_0 + c_i*tempU_1 + d_i*tempU_2)*tempFactor;
    //                         coefSB2[countB] = c_k*(b_i*tempU_0 + c_i*tempU_1 + d_i*tempU_2)*tempFactor;
    //                         coefSB3[countB] = d_k*(b_i*tempU_0 + c_i*tempU_1 + d_i*tempU_2)*tempFactor;
    //                         i_sparse[countB] = iglob;
    //                         k_sparse[countB] = kglob;
    //                         countB++;
    //                     }
    //                 }
    //             }
    //         }
    //         //----------------------------------
    //         // ADD SB1 and SB2 to SYSMAT_final          

    //         for (int locPos = 0; locPos < countB; locPos++)
    //         {
    //             int iglob = i_sparse[locPos];
    //             int kglob = k_sparse[locPos] + dim*nNodes_v;
    //             int pos = SYSMAT_final.getPos(iglob, kglob);
    //             coef[pos] += coefSB1[pos];
    //             coef[pos+rescaleCoef] += coefSB2[pos];
    //             coef[pos+2*rescaleCoef] += coefSB3[pos];
    //         }
    //         break;
    //     }
    // }
    // APPLY CORRECTION TO THE FORCING TERM
    switch (dim)
    {
        case 2:
        {
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* tempElem = elem_v[iel];
                std::vector<VECTOR> beta(dim); // each  beta[i] has the VECTOR with the i_th component of velocity in the nodes
                VECTOR f_sum(dim);


                prec tempFactor = Volume_v[iel]/12 * tau[iel];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    beta[icomp].initialize(dim+1);
                    std::shared_ptr<prec[]> tempP = beta[icomp].P;

                    prec* U_r_comp = U_r[icomp];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = tempElem[iloc];
                        tempP[iloc] = U_r_comp[iglob];

                        int globId = iglob+icomp*nNodes_v;
                        f_sum[icomp] += currForcingAtNodes[globId];
                    }
                }
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    prec b_i = Bloc_v[iel][iloc];
                    prec c_i = Cloc_v[iel][iloc];
                    int iglob = tempElem[iloc];

                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        int globId = iglob+icomp*nNodes_v;
                        prec tempF = currForcingAtNodes[globId] + f_sum[icomp];


                        rhs_final[globId] += (b_i*beta[0][iloc]*tempF + c_i*beta[1][iloc]*tempF) * tempFactor;
                    }
                }
            }
            break;
        }
        case 3:
        {
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* tempElem = elem_v[iel];
                std::vector<VECTOR> beta(dim); // each  beta[i] has the VECTOR with the i_th component of velocity in the nodes
                VECTOR f_sum(dim);


                prec tempFactor = Volume_v[iel]/20 * tau[iel];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    beta[icomp].initialize(dim+1);
                    std::shared_ptr<prec[]> tempP = beta[icomp].P;

                    prec* U_r_comp = U_r[icomp];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = tempElem[iloc];
                        tempP[iloc] = U_r_comp[iglob];

                        int globId = iglob+icomp*nNodes_v;
                        f_sum[icomp] += currForcingAtNodes[globId];
                    }
                }
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    prec b_i = Bloc_v[iel][iloc];
                    prec c_i = Cloc_v[iel][iloc];
                    prec d_i = Dloc_v[iel][iloc];
                    int iglob = tempElem[iloc];

                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        int globId = iglob+icomp*nNodes_v;
                        prec tempF = currForcingAtNodes[globId] + f_sum[icomp];


                        rhs_final[globId] += (b_i*beta[0][iloc] + c_i*beta[1][iloc] + d_i*beta[2][iloc]) * tempF * tempFactor;
                    }
                }
            }
            break;
        }
    }
}

//---------------------------------------
// BUILD Ma OPTIMIZATION PARAMETER MATRIX
//---------------------------------------
void PROBLEM_NS::assembleMa()
{
    // parameters copy;
    int dim = (*physics).dim;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;

    int nVar = dim+1;
    int tempQuotientBase = 10*(dim-1);
    int count = 0;
    VECTOR coefMa(nEvals);
    // assemble
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int* tempEl = elem_v[iel];
        prec tempVol = Volume_v[iel];
        for (int iloc = 0; iloc < nVar; iloc++)
        {
            // int iglob = tempEl[iloc];
            for (int jloc = 0; jloc < nVar; jloc++)
            {
                // int jglob = tempEl[jloc];
                prec tempVal = 0;
                int deltaOld = 0;
                prec tempQuotient = tempQuotientBase;
                if (iloc == jloc) deltaOld++;
                for (int kloc = 0; kloc < nVar; kloc++)
                {
                    int kglob = tempEl[kloc];
                    int delta = deltaOld;
                    if (kloc == iloc) delta++;
                    if (kloc == jloc) delta++;
                    switch (delta)
                    {
                        case 0:
                        {
                            tempQuotient *= 6;
                            break;
                        }
                        case 1:
                        {
                            tempQuotient *= 3;
                            break;
                        }
                    }
                    tempVal += alpha[kglob] * tempVol / tempQuotient;
                }
                coefMa[count] = tempVal;
                count++;
            }
        }
    }
    // resize coef to match H pattern
    resizeCoef(coefMa.P, Ma);
    // dlt
}

//--------------------------------------
// RESIZE COEF
//--------------------------------------
void PROBLEM_NS::resizeCoef(std::shared_ptr<prec[]> &coef, VECTOR &resized)
{
    resized.resetZeros();
    for (int ieval = 0; ieval < nEvals; ieval++)
    {
        resized[realPos[ieval]] += coef[ieval];
    }
}

//------------------------------------------------------
// ADD MATRIX TO H PATTERN IN SYSMAT BASE FOR NS PROBLEM
//------------------------------------------------------
void PROBLEM_NS::addToSysmat(std::shared_ptr<prec[]> &coef, prec factor)
{
    addToSysmat(coef, SYSMAT, factor);
}
//
void PROBLEM_NS::addToSysmat(std::shared_ptr<prec[]> &coef, CSRMAT &mat, prec factor)
{
    int dim = (*physics).dim;
    std::shared_ptr<prec[]> sysmatCoef = mat.coef;
    if (abs(factor - 1) < 1e-15)
    {
        #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int iterm = 0; iterm < nTerms; iterm++)
            {       
                int tempPos = realPosSysmat[0][icomp][iterm];
                sysmatCoef[tempPos] += coef[iterm]; 
            }
        }
    }
    else
    {
        #pragma omp parallel for num_threads(PARALLEL::nThread)
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int iterm = 0; iterm < nTerms; iterm++)
            {       
                int tempPos = realPosSysmat[0][icomp][iterm];
                sysmatCoef[tempPos] += factor * coef[iterm]; 
            }
        }
    }   
}
//---

