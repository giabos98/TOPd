#include "CODE_HEADERS/codeHeader.h"
#include "adjointNS.h"

void ADJOINT_NS::assembleA(VECTOR &velocity)
{
    // initialize geometry entities to not spend timme accessing to object parameters every time
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    //---

    // useful parameters
    int nLocNod = dim+1;
    int count;
    //---

    // create vector of shape functions coefs in prevision to have an already created one
    std::vector<MATRIX> coefMat(dim);
    switch (dim)
    {
        case 2:
        {
            coefMat[0] = (*physics).Bloc_v;
            coefMat[1] = (*physics).Cloc_v;
            break;
        }
        case 3:
        {
            coefMat[0] = (*physics).Bloc_v;
            coefMat[1] = (*physics).Cloc_v;
            coefMat[2] = (*physics).Dloc_v;
            break;
        }  
    }
    //---

    // build velocity components
    std::vector<VECTOR> vel(dim);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        vel[icomp].initialize(nNodes_v);
        for (int inod = 0; inod < nNodes_v; inod++)
        {
            vel[icomp][inod] = velocity[icomp*nNodes_v + inod];
        }
    }

    //------------------------------------------------------------------------
    // BUILD Ac u.dot(grad(u_a) over Volume & u.dot(n)ua over Neumann Boundary
    //------------------------------------------------------------------------
    VECTOR coefAc(nEvals);
    // coefAc = (prec*) malloc(nEvals * sizeof(prec)); checkPtr(coefAc);
    count = 0;
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int* globEl = elem_v[iel];
        prec tempFactor = Volume_v[iel]/nLocNod;

        // build mean velocity over elements per edge
        std::vector<VECTOR> locVel(dim);
        std::vector<VECTOR> locMeanVel(dim);
        VECTOR onesToSum(nLocNod); // vector to get the sum of a vector components doing scalar product
        onesToSum.reset(1.0);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            locVel[icomp].initialize(nLocNod);
            locMeanVel[icomp].setZeros(nLocNod);
            for (int inod = 0; inod < nLocNod; inod++)
            {
                int iglob = globEl[inod];
                prec tempVel = vel[icomp][iglob];
                locVel[icomp][inod] = tempVel;
                locMeanVel[icomp][inod] += tempVel;
            }
            prec tempSum = locVel[icomp].dot(onesToSum);
            locMeanVel[icomp] += tempSum;
            locMeanVel[icomp] /= (nLocNod+1);
        }
        //---

        // assemble matrix
        for (int iloc = 0; iloc < nLocNod; iloc++)
        {
            //int iglob = globEl[iloc];
            for (int jloc = 0; jloc < nLocNod; jloc++)
            {
                //int jglob = globEl[jloc];
                VECTOR tempMean(dim);
                VECTOR tempLocalBasis(dim);
                for (int idim = 0; idim < dim; idim++)
                {
                    // integrating -u.dot(grad(ua))) over the Volume and the Neumann part: u.dot(n)ua toghether
                    // in order to get implicit Neumann BC also in the adjoint weak formulation 
                    tempMean[idim] = locMeanVel[idim][jloc];  
                    tempLocalBasis[idim] = coefMat[idim][iel][iloc];
                } 
                coefAc[count] = tempMean.dot(tempLocalBasis) * tempFactor;
                count++;
            }
        }
        //---
    }
    std::shared_ptr<prec[]> AcP = coefAc.P;
    resizeCoef(AcP, Ac);

    //---

    //---
    
    //-----------------------------
    // BUILD Ag // grad(u).dot(u_a)
    //-----------------------------
    VECTOR coefAg(nEvals);
    prec tempQuotient = (dim+1)*(dim+2);
    for (int kdim = 0; kdim < dim; kdim++) // rows: determined by the test function component
    {                    
        for (int kcomp = 0; kcomp < dim; kcomp++) // cols: determined by the ua component
        {
            count = 0;
            coefAg.resetZeros();            
            for (int iel = 0; iel < nElem_v; iel++)
            {
                int* globEl = elem_v[iel];
                int tempFactor = Volume_v[iel] / tempQuotient;

                //create local velocity components
                VECTOR locVel(nLocNod);
                for (int inod = 0; inod < nLocNod; inod++)
                {
                    locVel[inod] = vel[kdim][globEl[inod]];
                }           
                //---

                // assemble
                for (int iloc = 0; iloc < nLocNod; iloc++)
                {
                    //int iglob = globEl[iloc];
                    for (int jloc = 0; jloc < nLocNod; jloc++)
                    {
                        //int jglob = globEl[jloc];
                        // build
                        prec tempFactor_ij = tempFactor;
                        if (iloc == jloc) tempFactor_ij *= 2;
                        prec tempAg = 0;
                        for (int knod = 0; knod < nLocNod; knod++)
                        { 
                            tempAg += locVel[knod]*coefMat[kdim][iel][knod] * tempFactor;
                        }
                        // coefAg[kdim][kcomp][count] = tempAg[kcomp];    
                        coefAg[count] = tempAg;                 
                        count++;                  
                    } 
                }
                //---
            }  
            resizeCoef(coefAg.P, Ag[kdim][kcomp]);
            //---
        }
    }
    //---
}
//---

//-----------------------------------------------------------------------------
// SYSMAT_AJD BASE
//-----------------------------------------------------------------------------
void ADJOINT_NS::createSysmatBase(CSRMAT &H, std::vector<CSRMAT> &B)
{
    // initialize adjoint matrices
    nTerms = H.nTerm;
    int dim = (*physics).dim;
    Ac.initialize(nTerms);
    
    Ag.resize(dim);
    for (int idim = 0; idim < dim; idim++)
    {
        Ag[idim].resize(dim);
        for (int jdim = 0; jdim < dim; jdim++)
        {
            Ag[idim][jdim].initialize(nTerms);
        }
    }

    Ma.initialize(nTerms);
    //----------------------
    int nNodes = (*physics).nNodes;
    int nNodes_v = (*physics).nNodes_v;

    int nTerm = dim * (dim*H.nTerm + 2*B[0].nTerm);
    int BLastRow = nNodes;
    if (nNeuBound == 0)
    {
        BLastRow = nNodes-1;
    }

    realPosSysmat.resize(dim);
    for (int idim = 0; idim < dim; idim++)
    {
         realPosSysmat[idim].resize(dim);
        for (int jdim = 0; jdim < dim; jdim++)
        {
            realPosSysmat[idim][jdim].initialize(nTerms);
        }
    } 

    int n = nNodes_v*dim + nNodes;
    SYSMAT_base.initialize(n, n, nTerm);
    std::shared_ptr<int[]> iat    = SYSMAT_base.iat;
    iat[0] = 0;
    std::shared_ptr<int[]> ja     = SYSMAT_base.ja;
    std::shared_ptr<prec[]> coef  = SYSMAT_base.coef;

    std::shared_ptr<int[]> iatH  = H.iat;
    std::shared_ptr<int[]>jaH   = H.ja;
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
            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                for (int ipos_v = iatH[irow]; ipos_v < iatH[irow+1]; ipos_v++)
                {
                    prec tempPos = jaH[ipos_v] + jcomp*nNodes_v;
                    ja[pos] = tempPos;
                    if (icomp == jcomp) coef[pos] = coefH[ipos_v];
                    else coef[pos] = 0;
                    realPosSysmat[icomp][jcomp][ipos_v] = pos;
                    pos++;
                }
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
        // for (int inode = 0; inode < nNodes; inode++)
        // {
        //     ja[pos]   = dim*nNodes_v + inode;
        //     coef[pos] = P_bc[inode];
        //     pos++;
        // }
        globRow++;
        iat[globRow] = pos;
    }

    if (nSymm > 0) 
    {
        CSRMAT tempMat = SYSMAT_base + J;
        SYSMAT_base = tempMat;
    }
}

//--------------------------------------------------------------
// ADD MATRIX TO H PATTERN IN SYSMAT BASE OF THE ADJOINT PROBLEM
//--------------------------------------------------------------
void ADJOINT_NS::addToSysmat(std::vector<std::vector<VECTOR>> &coef, int selectPath, prec factor)
{
    int dim = (*physics).dim;
    std::shared_ptr<prec[]> sysmatCoef = SYSMAT.coef;

    switch (selectPath)
    {
    case 0: // add only to diagonal blocks
    {
        std::shared_ptr<prec[]> coefPtr = coef[0][0].P;
        if ( abs(factor - 1) < 1e-15)
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int iterm = 0; iterm < nTerms; iterm++)
                {       
                    int tempPos = realPosSysmat[icomp][icomp][iterm];
                    sysmatCoef[tempPos] += coefPtr[iterm];
                }
            }
        }
        else
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int iterm = 0; iterm < nTerms; iterm++)
                {       
                    int tempPos = realPosSysmat[icomp][icomp][iterm];
                    sysmatCoef[tempPos] += factor * coefPtr[iterm]; 
                }
            }
        } 
        break;
    }
    case 1: // add to all blocks
    {
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                std::shared_ptr<prec[]> coefPtr = coef[icomp][jcomp].P;
                if ( abs(factor - 1) < 1e-15)
                {
                    for (int iterm = 0; iterm < nTerms; iterm++)
                    {       
                        int tempPos = realPosSysmat[icomp][jcomp][iterm];
                        sysmatCoef[tempPos] += coefPtr[iterm];
                    }
                }
                else
                {
                    for (int iterm = 0; iterm < nTerms; iterm++)
                    {       
                        int tempPos = realPosSysmat[icomp][jcomp][iterm];
                        sysmatCoef[tempPos] += factor * coefPtr[iterm]; 
                    }
                }
            }
        }
        break;
    }   
    }  
}
//---

void ADJOINT_NS::assembledAdgu()
{
    // int dim = (*physics).dim;
    // prec ni = (*physics).ni;
    // int nElem_v = (*physics).nElem_v;
    // std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    // std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;

    // std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    // std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    // std::shared_ptr<prec*[]> Dloc_v = (*physics).Dloc_v.PP;

    // int nNodxEl = dim+1;
    // int count = 0;
    // int nTermH = H.nTerm;
    
    // std::vector<std::vector<VECTOR*>> adjRHS(dim);
    // for (int i = 0; i < dim; i++) 
    // {
    //     adjRHS[i].resize(dim);
    // }
    // VECTOR coefH; coefH.length = nTermH; coefH.P = H.coef;

    // VECTOR mat00(nTermH); 
    // VECTOR mat01(nTermH);  
    // VECTOR mat02(nTermH); 
    // VECTOR mat10(nTermH); 
    // VECTOR mat11(nTermH);  
    // VECTOR mat12(nTermH); 
    // VECTOR mat20(nTermH); 
    // VECTOR mat21(nTermH);  
    // VECTOR mat22(nTermH); 

    // switch (dim)
    // {
    //     case 2:
    //     {
    //         //-----
            
    //         std::shared_ptr<prec[]> coef00(new prec[nEvals]);   
    //         std::shared_ptr<prec[]> coef01(new prec[nEvals]);  
    //         std::shared_ptr<prec[]> coef10(new prec[nEvals]);   
    //         std::shared_ptr<prec[]> coef11(new prec[nEvals]);  
    //         for (int iel = 0; iel < nElem_v; iel++)
    //         {

    //             prec tempFactor = ni*Volume_v[iel];
    //             for (int iloc = 0; iloc < nNodxEl; iloc++)
    //             {
    //                 prec b_i = Bloc_v[iel][iloc];
    //                 prec c_i = Cloc_v[iel][iloc];
    //                 for (int jloc = 0; jloc < nNodxEl; jloc++)
    //                 {
    //                     prec b_j = Bloc_v[iel][jloc];
    //                     prec c_j = Cloc_v[iel][jloc];
    //                     prec fact = 0;
    //                     if (customFunc == 0)
    //                     {
    //                         if (onlyGrad == 0) fact = 0; 
    //                         if (iel == 0 && iloc == 0 && jloc == 0) std::cout << "\nqua\n";
    //                         coef00[count] = fact*(b_j*b_i) * tempFactor;
    //                         coef01[count] = fact*(b_j*c_i) * tempFactor;
    //                         coef10[count] = fact*(c_j*b_i) * tempFactor;
    //                         coef11[count] = fact*(c_j*c_i) * tempFactor;

    //                     }
    //                     else if (customFunc == 1)
    //                     {
    //                         // gradT
    //                         if (onlyGrad == 0) fact = fWeights[1]; 

    //                         coef00[count] = fact*(b_j*b_i) * tempFactor;
    //                         coef01[count] = fact*(b_j*c_i) * tempFactor;
    //                         coef10[count] = fact*(c_j*b_i) * tempFactor;
    //                         coef11[count] = fact*(c_j*c_i) * tempFactor;
                            
    //                         // curl
    //                         coef00[count] += fWeights[2]*(c_j*c_i) * tempFactor;
    //                         coef01[count] += fWeights[2]*(-b_j*c_i) * tempFactor;
    //                         coef10[count] += fWeights[2]*(-c_j*b_i) * tempFactor;
    //                         coef11[count] += fWeights[2]*(b_j*b_i) * tempFactor;
    //                     }   
    //                     count ++;
    //                 }
    //             }
    //         }
    //         resizeCoef(coef00, mat00); 
    //         resizeCoef(coef01, mat01); 
    //         resizeCoef(coef10, mat10);
    //         resizeCoef(coef11, mat11);  

    //         if (customFunc == 0)
    //         {
    //             PARALLEL::sum(mat00, coefH, mat00);
    //             PARALLEL::sum(mat11, coefH, mat11);
    //         }
    //         else if (customFunc == 1)
    //         {
    //             VECTOR coefHfinal(nTermH); coefH.multiplyByScalar(fWeights[1], coefHfinal);
    //             PARALLEL::sum(mat00, coefHfinal, mat00);
    //             PARALLEL::sum(mat11, coefHfinal, mat11);
    //         }
            

    //         adjRHS[0][0] = &mat00;  adjRHS[0][1] = &mat01;
    //         adjRHS[1][0] = &mat10;  adjRHS[1][1] = &mat11;    
    //         break;
    //     }
    //     case 3:
    //     {
    //         //-----
            
    //         std::shared_ptr<prec[]> coef00(new prec[nEvals]);   
    //         std::shared_ptr<prec[]> coef01(new prec[nEvals]); 
    //         std::shared_ptr<prec[]> coef02(new prec[nEvals]);  
    //         std::shared_ptr<prec[]> coef10(new prec[nEvals]);   
    //         std::shared_ptr<prec[]> coef11(new prec[nEvals]);  
    //         std::shared_ptr<prec[]> coef12(new prec[nEvals]);  
    //         std::shared_ptr<prec[]> coef20(new prec[nEvals]); 
    //         std::shared_ptr<prec[]> coef21(new prec[nEvals]); 
    //         std::shared_ptr<prec[]> coef22(new prec[nEvals]);   
    //         for (int iel = 0; iel < nElem_v; iel++)
    //         {

    //             prec tempFactor = ni*Volume_v[iel];
    //             for (int iloc = 0; iloc < nNodxEl; iloc++)
    //             {
    //                 prec b_i = Bloc_v[iel][iloc];
    //                 prec c_i = Cloc_v[iel][iloc];
    //                 prec d_i = Dloc_v[iel][iloc];
    //                 for (int jloc = 0; jloc < nNodxEl; jloc++)
    //                 {
    //                     prec b_j = Bloc_v[iel][jloc];
    //                     prec c_j = Cloc_v[iel][jloc];
    //                     prec d_j = Dloc_v[iel][jloc];
    //                     prec fact = 0;
    //                     if (customFunc == 0)
    //                     {
    //                         if (onlyGrad == 0) fact = 0; 
    //                         coef00[count] = fact*(b_j*b_i) * tempFactor;
    //                         coef01[count] = fact*(b_j*c_i) * tempFactor;
    //                         coef02[count] = fact*(b_j*d_i) * tempFactor;
    //                         coef10[count] = fact*(c_j*b_i) * tempFactor;
    //                         coef11[count] = fact*(c_j*c_i) * tempFactor;
    //                         coef12[count] = fact*(c_j*d_i) * tempFactor;
    //                         coef20[count] = fact*(d_j*b_i) * tempFactor;
    //                         coef21[count] = fact*(d_j*c_i) * tempFactor;
    //                         coef22[count] = fact*(d_j*d_i) * tempFactor;
    //                     }
    //                     else if (customFunc == 1)
    //                     {
    //                         // gradT
    //                         if (onlyGrad == 0) fact = fWeights[1]; 

    //                         coef00[count] = fact*(b_j*b_i) * tempFactor;
    //                         coef01[count] = fact*(b_j*c_i) * tempFactor;
    //                         coef02[count] = fact*(b_j*d_i) * tempFactor;
    //                         coef10[count] = fact*(c_j*b_i) * tempFactor;
    //                         coef11[count] = fact*(c_j*c_i) * tempFactor;
    //                         coef12[count] = fact*(c_j*d_i) * tempFactor;
    //                         coef20[count] = fact*(d_j*b_i) * tempFactor;
    //                         coef21[count] = fact*(d_j*c_i) * tempFactor;
    //                         coef22[count] = fact*(d_j*d_i) * tempFactor;
                            
    //                         // curl
    //                         coef00[count] += fWeights[2]*(c_j*c_i + d_j*d_i) * tempFactor;
    //                         coef01[count] += fWeights[2]*(-b_j*c_i) * tempFactor;
    //                         coef02[count] += fWeights[2]*(-b_j*d_i) * tempFactor;
    //                         coef10[count] += fWeights[2]*(-c_j*b_i) * tempFactor;
    //                         coef11[count] += fWeights[2]*(d_j*d_i + b_j*b_i) * tempFactor;
    //                         coef12[count] += fWeights[2]*(-c_j*d_i) * tempFactor;
    //                         coef20[count] += fWeights[2]*(-d_j*b_i) * tempFactor;
    //                         coef21[count] += fWeights[2]*(-d_j*c_i) * tempFactor;
    //                         coef22[count] += fWeights[2]*(b_j*b_i + c_j*c_i) * tempFactor;
    //                     }   
    //                     count ++;
    //                 }
    //             }
    //         }
    //         resizeCoef(coef00, mat00); 
    //         resizeCoef(coef01, mat01); 
    //         resizeCoef(coef02, mat02); 
    //         resizeCoef(coef10, mat10);
    //         resizeCoef(coef11, mat11); 
    //         resizeCoef(coef12, mat12); 
    //         resizeCoef(coef20, mat20);            
    //         resizeCoef(coef21, mat21);
    //         resizeCoef(coef22, mat22); 

    //         if (customFunc == 0)
    //         {
    //             PARALLEL::sum(mat00, coefH, mat00);
    //             PARALLEL::sum(mat11, coefH, mat11);
    //             PARALLEL::sum(mat22, coefH, mat22);
    //         }
    //         else if (customFunc == 1)
    //         {
    //             VECTOR coefHfinal(nTermH); coefH.multiplyByScalar(fWeights[1], coefHfinal);
    //             PARALLEL::sum(mat00, coefHfinal, mat00);
    //             PARALLEL::sum(mat11, coefHfinal, mat11);
    //             PARALLEL::sum(mat22, coefHfinal, mat22);
    //         }
            

    //         adjRHS[0][0] = &mat00;  adjRHS[0][1] = &mat01;  adjRHS[0][2] = &mat02;  
    //         adjRHS[1][0] = &mat10;  adjRHS[1][1] = &mat11;  adjRHS[1][2] = &mat12;  
    //         adjRHS[2][0] = &mat20;  adjRHS[2][1] = &mat21;  adjRHS[2][2] = &mat22;  
    //         break;
    //     }
    // }

    // PARALLEL::createBlockMatrix(H, adjRHS, dAdGu);

    // prec two = -2.0;
    // std::shared_ptr<prec[]> coef = dAdGu.coef;
    // VECTOR::multiplyByScalar(coef, two, dAdGu.nTerm, coef);

    int dim = (*physics).dim;
    prec ni = (*physics).ni;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;

    std::vector<MATRIX*> Coef_v(dim);
    for (int icomp = 0 ; icomp < dim; icomp++)
    {
        Coef_v[icomp] = &(*physics).Coef_v[icomp];
    }

    int nNodxEl = dim+1;
    int count = 0;
    int nTermH = H.nTerm;

    prec weight_grad   = -(2*fWeights[1] + fWeights[2]);
    prec weigth_grad_t = -(2*fWeights[1] - fWeights[2]);
    

    VECTOR coefH; coefH.length = nTermH; coefH.P = H.coef;

    std::vector<std::vector<VECTOR>> mat_bloc_vect(dim);
    std::vector<std::vector<VECTOR>> mat_bloc_vect_coef(dim);
    std::vector<std::vector<VECTOR*>> adjRHS(dim);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        mat_bloc_vect[icomp].resize(dim);
        mat_bloc_vect_coef[icomp].resize(dim);
        adjRHS[icomp].resize(dim);
        for (int jcomp = 0; jcomp < dim; jcomp++)
        {
            mat_bloc_vect[icomp][jcomp].initialize(nTermH);  
            mat_bloc_vect[icomp][jcomp].resetZeros();  
            mat_bloc_vect_coef[icomp][jcomp].initialize(nEvals);  
            mat_bloc_vect_coef[icomp][jcomp].resetZeros();  
        }
    }
 
    //pause();

    for (int iel = 0; iel < nElem_v; iel++)
    {
        prec tempFactor = ni*Volume_v[iel];
        for (int iloc = 0; iloc < nNodxEl; iloc++)
        {
            for (int jloc = 0; jloc < nNodxEl; jloc++)
            {
                //std::cout << "iel: " << iel << "\tiloc: " << iloc << "\tjloc: " << jloc << "\n";
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    //std::cout << "\t\ticomp: " << icomp << "\n";
                    //prec coef_icomp_iloc = *(Coef_v[icomp])[iel][iloc];
                    prec coef_icomp_jloc = (*physics).Coef_v[icomp][iel][jloc];
                    for (int jcomp = 0; jcomp < dim; jcomp++)
                    {
                        //std::cout << "\t\t\tjcomp: " << jcomp << "\n";
                        prec coef_jcomp_iloc = (*physics).Coef_v[icomp][iel][iloc];
                        //prec coef_jcomp_jloc = *(Coef_v[jcomp])[iel][jloc];
                        //std::cout << "\t\t\t\t ASSIGN\n";
                        mat_bloc_vect_coef[icomp][jcomp][count] += weigth_grad_t * (coef_icomp_jloc * coef_jcomp_iloc) * tempFactor;
                        //std::cout << "\t\t\t\t DONE\n";
                    }
                }
                count ++;
            }
        }
    }
    //pause();

    for (int icomp = 0; icomp < dim; icomp++)
    {
        for (int jcomp = 0; jcomp < dim; jcomp++)
        {
            resizeCoef(mat_bloc_vect_coef[icomp][jcomp].P, mat_bloc_vect[icomp][jcomp]);
        }
    }

    VECTOR coefHfinal(nTermH); coefH.multiplyByScalar(weight_grad, coefHfinal);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        PARALLEL::sum(mat_bloc_vect[icomp][icomp], coefHfinal, mat_bloc_vect[icomp][icomp]);
    }
    
    for (int icomp = 0; icomp < dim; icomp++)
    {
        for (int jcomp = 0; jcomp < dim; jcomp++)
        {
            adjRHS[icomp][jcomp] = &mat_bloc_vect[icomp][jcomp];
        }
    }        

    PARALLEL::createBlockMatrix(H, adjRHS, dAdGu);
    //pause();
}