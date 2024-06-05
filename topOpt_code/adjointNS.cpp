#include "adjointNS.h"

//----------------------------------
// INITIALIZE
//----------------------------------
 void ADJOINT_NS::initialize(PROBLEM_NS &NS, VECTOR &alphaIn)
 {
    physics = NS.physics;
    name = NS.name;
    time = 0; 
    deltaT = NS.deltaT;
    save_solution = NS.save_solution;
    completeLog = NS.completeLog;

    alpha.length = alphaIn.length;
    alpha.P = alphaIn.P;
    nDof = NS.nDof;

    std::string folderName;
    switch ((*physics).isNS)
    {
        case 0:
        {
            folderName = "ADJ_S_sol";
            break;
        }
        case 1:
        {
            folderName = "ADJ_NS_sol";
            break;
        }
        default:
        {
            break;
        }
    }
    VTKWriter.initializeForADJ(name, (*physics).dim, (*physics).nNodes, (*physics).nElem, folderName);
}
//----------------------------------
// SET ADJOINT BC
//-----------------------------------
void ADJOINT_NS::setBC(PROBLEM_NS &NS)
{
    printf("\n-------\n-| SET ADJ BC |--\n-------\n");
    //---
    // access NS parameters
    //- Wall BC --- (priority: 5)
    int nWallBoundNS = NS.nWallBound;
    VECTOR_INT wallBoundNS;
    wallBoundNS = NS.wallBound;
    VECTOR_INT wallNodNS;
    wallNodNS = NS.wallNod; 
    VECTOR_INT wallIdCountNS;
    wallIdCountNS = NS.wallIdCount;
    int  nWallNS = NS.nWall;

    //- Temporal Dirichlet BC --- (priority: 4)
    int nDirTimeBoundNS = NS.nDirTimeBound;
    VECTOR_INT dirTimeBoundNS;
    dirTimeBoundNS = NS.dirTimeBound;
    VECTOR_INT dirTimeNodNS;
    dirTimeNodNS = NS.dirTimeNod; 
    VECTOR_INT dirTimeIdCountNS;
    dirTimeIdCountNS = NS.dirTimeIdCount;
    int nTimeDirNS = NS.nTimeDir;

    //- Static Dirichlet BC --- (priority: 3)
    int nDirBoundNS = NS.nDirBound;
    VECTOR_INT dirBoundNS;
    dirBoundNS = NS.dirBound;
    VECTOR_INT dirNodNS;
    dirNodNS = NS.dirNod; 
    VECTOR_INT dirIdCountNS;
    dirIdCountNS = NS.dirIdCount;
    
    int  nDirNS = NS.nDir;
    
    //- Temporal Neumann BC --- (priority: 2)
    int nNeuTimeBoundNS = NS.nNeuTimeBound;
    VECTOR_INT neuTimeBoundNS;
    neuTimeBoundNS = NS.neuTimeBound;
    VECTOR_INT neuTimeNodNS;
    neuTimeNodNS = NS.neuTimeNod; 
    VECTOR_INT neuTimeIdCountNS;
    neuTimeIdCountNS = NS.neuTimeIdCount;
    int nTimeNeuNS = NS.nTimeNeu;
    
    //- Static Neumann BC --- (priority: 1)
    int nNeuBoundNS = NS.nNeuBound;
    VECTOR_INT neuBoundNS;
    neuBoundNS = NS.neuBound;
    VECTOR_INT neuNodNS;
    neuNodNS = NS.neuNod;
    VECTOR_INT neuIdCountNS;
    neuIdCountNS = NS.neuIdCount;
    int nNeuNS = NS.nNeu;

    //- NVOB BC --- (priority: -2)
    VECTOR_INT NVOBBoundNS;
    NVOBBoundNS = NS.NVOBBound;
    VECTOR_INT NVOBNodNS;
    NVOBNodNS = NS.NVOBNod;
    VECTOR_INT NVOBIdCountNS;
    NVOBIdCountNS = NS.NVOBIdCount;

    int NodesxBoundElem = (*physics).dim;
    bool isStat = 1; // parameter to describe is the problem is stationary or not
    // Dirichlet ADJ BC
    if (abs(fWeights[3]) < 1e-12) // dir -> wall
    {
        nWall      = nWallNS + nTimeDirNS + nDirNS;
        nWallBound = nWallBoundNS + nDirTimeBoundNS + nDirBoundNS;
        wallBound  = wallBoundNS;
        wallBound.append(dirTimeBoundNS);
        wallBound.append(dirBoundNS);
        //---
        wallIdCount.initialize(nWallBound);
        int iwallBound = 0;
        for (int iwallNS = 0; iwallNS < nWallBoundNS; iwallNS++)
        {
            wallIdCount[iwallBound] = wallIdCountNS[iwallNS]; iwallBound++;
        }
        for (int idirTimeNS = 0; idirTimeNS < nDirTimeBoundNS; idirTimeNS++)
        {
            wallIdCount[iwallBound] = dirTimeIdCountNS[idirTimeNS] + nWallNS; iwallBound++;
        } 
        for (int idirNS = 0; idirNS < nDirBoundNS; idirNS++)
        {
            wallIdCount[iwallBound] = dirIdCountNS[idirNS] + nWallNS + nTimeDirNS; iwallBound++;
        } 
        //---
        wallNod = wallNodNS;
        wallNod.append(dirTimeNodNS);
        wallNod.append(dirNodNS);
        //---
        nDir = 0;
        nDirBound = 0;
        nTimeDir = 0;
        nDirTimeBound = 0;
    }
    else // wall -> dir
    {
        //if (isStat) // stationary COMMENTED CAUSE IN THE PRESSURE FUNCTIONAL, IN ANY CASE THE dB/dp = 1 constantly
        //{
        nDir          = nWallNS + nTimeDirNS + nDirNS;
        nDirBound     = nWallBoundNS + nDirTimeBoundNS + nDirBoundNS;
        dirBound      = wallBoundNS;
        dirBound.append(dirTimeBoundNS);
        dirBound.append(dirBoundNS);
        //---
        dirIdCount.initialize(nDirBound);
        int idirBound = 0;
        for (int iwallNS = 0; iwallNS < nWallBoundNS; iwallNS++)
        {
            dirIdCount[idirBound] = wallIdCountNS[iwallNS]; idirBound++;
        }
        for (int idirTimeNS = 0; idirTimeNS < nDirTimeBoundNS; idirTimeNS++)
        {
            dirIdCount[idirBound] = dirTimeIdCountNS[idirTimeNS] + nWallNS; idirBound++;
        } 
        for (int idirNS = 0; idirNS < nDirBoundNS; idirNS++)
        {
            dirIdCount[idirBound] = dirIdCountNS[idirNS] + nWallNS + nTimeDirNS; idirBound++;
        } 
        //---
        dirNod  = wallNodNS;
        dirNod.append(dirTimeNodNS);
        dirNod.append(dirNodNS);
        //---
        nWall = 0;
        nWallBound = 0;
        nTimeDir = 0;
        nDirTimeBound = 0;

        std::string folderPath    = name;
        folderPath = "PREPRO/PROBLEM_DATA/" + folderPath;
        std::string BoundElemFilePath_v = folderPath + "/BoundElems_v.txt";
        std::ifstream  BoundElemFile_v; BoundElemFile_v.open(BoundElemFilePath_v, std::ios::in | std::ios::binary); if (!BoundElemFile_v.is_open()) throw_line("ERROR, can't open input data file");
        //--------------------------------
        // READ BOUND INFO FROM BOUND FILE
        //--------------------------------
        std::string line;
        std::istringstream iss;
        std::getline(BoundElemFile_v, line);
        iss.str(line);
        iss >> nBoundElems; iss >> nBound;
        iss.clear();

        //----- alloc dirBoundIdElems ---
        normalBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
        normalBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
        normalBoundIdElems      = std::shared_ptr<int**[]> (new int**[nDirBound]);
        nNormalBoundIdElems     = std::shared_ptr<int[]> (new int[nDirBound]);
        //---
        int countN = 0; int countElNormal = 0;
        int countNormal = 0;
        //---
        for (int i = 0; i < nBound; i++)
        {
            int currID;
            int currNEl;
            std::getline(BoundElemFile_v, line);
            iss.str(line);
            iss >> currID; iss >> currNEl;
            iss.clear();
            if (dirBound.hasIn(currID))
            {
                normalBoundIdElems[countNormal] = &normalBoundIdElems_buff[countElNormal];
                for (int j = 0; j < currNEl; j++)
                {
                    int* p = &normalBoundIdElems_buff_buff[countN];
                    normalBoundIdElems_buff[countElNormal] = p;
                    STREAM::getRowVector(BoundElemFile_v, line, iss, p, NodesxBoundElem);
                    countN += NodesxBoundElem;
                    countElNormal++;
                }
                nNormalBoundIdElems[countNormal] = currNEl;
                countNormal++;
                std::getline(BoundElemFile_v, line);
            } 
            else
            {
                STREAM::getLines(BoundElemFile_v, line, currNEl+1);
                continue;
            }
        }
        BoundElemFile_v.close();
            // }
            // else // time dependent: never the case cause derivating the B functional by the pressure we have a constant
            // {
            //     nTimeDir          = nWallNS + nTimeDirNS + nDirNS;
            //     nDirTimeBound     = nWallBoundNS + nDirTimeBoundNS + nDirBoundNS;
                
            //     dirTimeBound  = wallBoundNS;
            //     dirTimeBound.append(dirTimeBoundNS);
            //     dirTimeBound.append(dirBoundNS);
            //     //---
            //     //---
            //     dirTimeIdCount.initialize(nDir);
            //     int idirBound = 0;
            //     for (int iwallNS = 0; iwallNS < nWallBoundNS; iwallNS++)
            //     {
            //         dirTimeIdCount[idirBound] = wallIdCountNS[iwallNS]; idirBound++;
            //     }
            //     for (int idirTimeNS = 0; idirTimeNS < nDirTimeBoundNS; idirTimeNS++)
            //     {
            //         dirTimeIdCount[idirBound] = dirTimeIdCountNS[idirTimeNS] + nWallNS; idirBound++;
            //     } 
            //     for (int idirNS = 0; idirNS < nDirBoundNS; idirNS++)
            //     {
            //         dirTimeIdCount[idirBound] = dirIdCountNS[idirNS] + nWallNS + nTimeDirNS; idirBound++;
            //     } 
            //     //---
            //     dirTimeNod = wallNodNS;
            //     dirTimeNod.append(dirTimeNodNS);
            //     dirTimeNod.append(dirNodNS);
            //     // dirTimeNod        = VECTOR_INT::merge(wallNodNS, VECTOR_INT::merge(dirTimeNodNS, dirNodNS));
            //     //---
            //     nWall = 0;
            //     nWallBound = 0;
            //     nDir = 0;
            //     nDirBound = 0;
            // }
    }

    // Neumann ADJ BC
    if (isStat) // stationary
    {
        nNeu       = nTimeNeuNS + nNeuNS; 
        nNeuBound  = nNeuTimeBoundNS + nNeuBoundNS;
        neuBound   = neuTimeBoundNS; neuBound.append(neuBoundNS);
        neuIdCount = neuTimeIdCountNS; neuIdCount.append(neuIdCountNS);
        neuNod     = neuTimeNodNS; neuNod.append(neuNodNS);
        //---
        nTimeNeu = 0;
        nNeuTimeBound = 0;
    }
    else // time dependent
    {
        nTimeNeu       = nTimeNeuNS + nNeuNS; 
        nNeuTimeBound  = nNeuTimeBoundNS + nNeuBoundNS;
        neuTimeBound   = neuTimeBoundNS; neuTimeBound.append(neuBoundNS);
        neuTimeIdCount = neuTimeIdCountNS; neuTimeIdCount.append(neuIdCountNS);
        neuTimeNod     = neuTimeNodNS; neuTimeNod.append(neuNodNS);
        //---
        nNeu = 0;
        nNeuBound = 0;
    }

    dirVal.initialize(nDir, (*physics).dim);
    dirVal.resetZeros();
    neuVal.zeros(nNeu, (*physics).dim);

    if (NS.nSymmBound > 0)
    {
        symmBound = NS.symmBound;
        symmNod = NS.symmNod;
        symmVal = NS.symmVal;
        nSymmBound = NS.nSymmBound;
        nSymm = NS.nSymm;
        symmRealIglob = NS.symmRealIglob;
        symmVal = NS.symmVal;

        symmPosToZero = NS.symmPosToZero;
        symmPosToCoef = NS.symmPosToCoef;

        J = NS.J;
    }


    // NVOB BC
    NVOBBound = NS.NVOBBound;
    NVOBNod = NS.NVOBNod;
    nNVOBBound = NS.nNVOBBound;
    nNVOB = NS.nNVOB;

    setStatDirBC();
    setStatNeuBC();
}
//-----------------------
//---
void ADJOINT_NS::setStatDirBC() //priority 3
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| SET STAT DIR BC  |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
    MATRIX matCoord(dim,dim); 
    std::vector<VECTOR> normals(nNodes_v);
    VECTOR_INT inlet_nodes_v((*physics).inlet_nodes_v.P, (*physics).inlet_nodes_v.length);
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        normals[inod].setZeros(dim);
    }
    // MATRIX normalsFull(nNodes_v,dim); normalsFull.resetZeros();
    // VECTOR nodalAreas(nNodes_v); nodalAreas.resetZeros();
    VECTOR tempNormal(dim);
    // prec rho = (*physics).rho;
    // get normals for Dir nodes
    for (int ibound = 0; ibound < nDirBound; ibound++)
    {
        //-------------------------------------
        for (int iel = 0; iel < nNormalBoundIdElems[ibound]; iel++)
        {
            int* el = normalBoundIdElems[ibound][iel]; 
            for (int inode = 0; inode < dim; inode++)
            {
                int iglob = el[inode];
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
            }
            //prec Area = getArea(matCoord, dim);
            getNormal(matCoord, tempNormal);
            for (int inod = 0; inod < dim; inod++)
            {
                int iglob = el[inod];
                // for (int icomp = 0; icomp < dim; icomp ++) 
                // {
                //     normalsFull[iglob][icomp] += tempNormal[icomp]*Area;
                //     nodalAreas[iglob] += Area;
                // }  
                normals[iglob] += tempNormal;
            } 
        }
    }
    for (int inod = 0; inod < nNodes_v; inod++)
    {
        prec temp_norm = normals[inod].norm();
        normals[inod] = normals[inod] / temp_norm;
    }
    //impose Dir BC equal to the normals
    for (int idir = 0; idir < nDir; idir++)
    {
        int iglob = dirNod[idir];
        // for (int icomp = 0; icomp < dim; icomp++)
        // {
        //     dirVal[idir][icomp] = -fWeights[3]*normalsFull[iglob][icomp]/nodalAreas[iglob];
        // }
        if (inlet_nodes_v.hasIn(iglob))
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                dirVal[idir][icomp] = -fWeights[3]*normals[iglob][icomp];
            }
        }
    }
    // pause();
    // MATRIX::printForMatlab(dirVal, "dir values");
    // pause();

}
//---
void ADJOINT_NS::setStatNeuBC() // priority 1
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| SET STAT NEU BC |--\n----------\n";
    int dim = (*physics).dim;
    // prec** coord_v = (*physics).coord_v.PP;
    // prec rho = (*physics).rho;
    MATRIX matCoord(dim,dim); 
    VECTOR normal(dim);

    // for (int ibound = 0; ibound < nNeuBound; ibound++)
    // {
    //     prec meanPByDim = neuMeanP[ibound]/dim;
    //     //-------------------------------------
    //     for (int iel = 0; iel < nNeuBoundIdElems[ibound]; iel++)
    //     {
    //         int* el = neuBoundIdElems[ibound][iel];
    //         for (int inode = 0; inode < dim; inode++)
    //         {
    //             int iglob = el[inode];
    //             for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
    //         }
    //         prec Area = getArea(matCoord, dim);
    //         getNormal(matCoord, normal);
    //         //---
    //         prec tempFactor = Area * meanPByDim;

    //         for (int iloc = 0; iloc < dim; iloc++)
    //         {
    //             int iglob = el[iloc];
    //             // std::cout << "iglob: " << iglob << "\n";
    //             if (priority[iglob] != 1) 
    //             {
    //                 continue;
    //             }
    //             int tempNeuNod = neuNodPos[iglob];
    //             VECTOR value = normal * tempFactor;
    //             prec* tempNeu = neuVal[tempNeuNod];
    //             for (int icomp = 0; icomp < dim; icomp++)
    //             {
    //                 tempNeu[icomp] += value[icomp]/rho;
    //             }
    //             value.dlt();
    //         }    
    //     }
    // }
    // matCoord.dlt();
}

//------------------------------------------------------------------------------
// IMPOSE BOUNDARY CONDITIONS
//------------------------------------------------------------------------------
void ADJOINT_NS::imposeBC(CSRMAT &SYSMAT_final, VECTOR &rhs)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| IMPOSE BC |--\n----------\n";
    imposeStaticBC(SYSMAT_final, rhs);
    imposeTimeBC(SYSMAT_final, rhs);
}

//------------------------------
// IMPOSE STATIC BC
//------------------------------
void ADJOINT_NS::imposeStaticBC(CSRMAT &SYSMAT_NS, VECTOR &rhs)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| IMPOSE STATIC BC  |--\n----------\n";
    int dim = (*physics).dim;
    int flagBC = (*physics).flagBC;
    int nNodes_v = (*physics).nNodes_v;
    // DIRICHLET BC
    switch (flagBC)
    {
    case 0:
    {
        prec penalty = 1e13;
        for (int icomp = 0; icomp < dim; icomp++)
        {
            // WALL dir
            for (int iwall = 0; iwall < nWall; iwall++)
            {
                int iglob = wallNod[iwall] + nNodes_v*icomp; 
                SYSMAT_NS(iglob,iglob, penalty);
                rhs[iglob] = 0;
            }
            // static dir
            for (int idir = 0; idir < nDir; idir++)
            {
                int iglob = dirNod[idir] + nNodes_v*icomp;  
                SYSMAT_NS(iglob,iglob) = penalty;
                rhs[iglob] = penalty * dirVal[idir][icomp];
            }
        }
        for (int icomp = 0; icomp < dim-1; icomp++)
        {
            // NVOB dir
            for (int iNVOB = 0; iNVOB < nNVOB; iNVOB++)
            {
                int iglob = NVOBNod[iNVOB] + nNodes_v*icomp; 
                SYSMAT_NS(iglob,iglob, penalty);
                rhs[iglob] = 0;
            }
        }
        // symm 
        for (int isymm = 0; isymm < nSymm; isymm++)
        {
            int iglob = symmNod[isymm]; 

            prec* tempValue = symmVal[isymm];

            int trueIglob = 0;

            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                if (tempValue[jcomp] != 0)
                {
                    trueIglob = iglob+jcomp*nNodes_v;
                    break;
                }
            }
            
            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                int jglob = iglob + jcomp*nNodes_v;
                SYSMAT_NS(trueIglob, jglob) = tempValue[jcomp];
            } 
            rhs[trueIglob] = 0;
        }
        break;
    }
    case 1:
    {
        
        if (!statBCAlreadyApplied)
        {
            // WALL
            // set to zero row & col entries
            int countPosToZero = 0; int countIglobToZero = 0; int countPosToOne = 0;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int iwall = 0; iwall < nWall; iwall++)
                {
                    int iglob = wallNod[iwall] + icomp*nNodes_v; 
                    VECTOR_INT rowInColToZero = SYSMAT_NS.setColInRowToZero(iglob, wallPosToZero, countPosToZero, wallPosToOne, countPosToOne);  
                    rhs[iglob] = 0;                                             
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int tempPos = SYSMAT_NS(rowInColToZero[i], iglob, 0);
                        wallPosToZero[countPosToZero] = tempPos; countPosToZero++;
                    }
                }
            }
            //---
            // NVOB
            // set to zero row & col entries
            countPosToZero = 0; countIglobToZero = 0; countPosToOne = 0;
            for (int icomp = 0; icomp < dim-1; icomp++)
            {
                for (int iNVOB = 0; iNVOB < nNVOB; iNVOB++)
                {
                    int iglob = NVOBNod[iNVOB] + icomp*nNodes_v; 
                    VECTOR_INT rowInColToZero = SYSMAT_NS.setColInRowToZero(iglob, NVOBPosToZero, countPosToZero, NVOBPosToOne, countPosToOne);  
                    rhs[iglob] = 0;                                             
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int tempPos = SYSMAT_NS(rowInColToZero[i], iglob, 0);
                        NVOBPosToZero[countPosToZero] = tempPos; countPosToZero++;
                    }
                }
            }
            //---
            // static dir
            // set to zero row & col entries
            countPosToZero = 0; countPosToOne = 0;
            std::shared_ptr<prec[]> coefMat = SYSMAT_NS.coef;
            dirCountToZero[0] = 0;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int idir = 0; idir < nDir; idir++)
                {
                    prec tempDirVal = dirVal[idir][icomp];
                    int tempDirNod  = dirNod[idir] + icomp*nNodes_v;
                    rhs[tempDirNod] = tempDirVal;
                    VECTOR_INT rowInColToZero = SYSMAT_NS.setColInRowToZero(tempDirNod, dirPosToZero, countPosToZero, dirPosToOne, countPosToOne);
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int iglob = rowInColToZero[i];
                        int pos = SYSMAT_NS.getPos(iglob, tempDirNod);
                        prec tempVal = coefMat[pos];
                        rhs[iglob] -= tempVal * tempDirVal;
                        coefMat[pos] = 0.0;
                        dirIglobToZero[countIglobToZero] = iglob; countIglobToZero++;
                        dirPosToZero[countPosToZero] = pos; countPosToZero++;

                        int tempDirCountId = icomp*nDir + idir + 1;
                        dirCountToZero[tempDirCountId] = countPosToZero;
                    }
                }
            }

            //---
            // SYMMETRIC BC
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            countPosToZero = 0; int countPosToCoef = 0;
            std::shared_ptr<int[]> iat = SYSMAT_NS.iat;
            std::shared_ptr<int[]> ja = SYSMAT_NS.ja;
            for (int isymm = 0; isymm < nSymm; isymm++)
            {
                prec* tempVal = symmVal[isymm];
                int iglob  = symmNod[isymm];

                int trueIglob = 0;
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    if (abs(tempVal[jcomp]) > 0.8)
                    {
                        trueIglob = iglob+jcomp*nNodes_v;
                        break;
                    }
                }

                // VECTOR::print(tempVal,3);
                symmRealIglob[isymm] = trueIglob;

                rhs[trueIglob] = 0;
                
                int startPos = iat[trueIglob]; int endPos = iat[trueIglob+1];

                int jglob = iglob; int jcomp = 0;

                for (int pos = startPos; pos < endPos; pos++)
                {
                    int tempCol = ja[pos];
                    if (tempCol == jglob)
                    {
                        coefMat[pos] = tempVal[jcomp];
                        symmPosToCoef[countPosToCoef] = pos;
                        countPosToCoef++;
                        if (jcomp != dim-1) 
                        {
                            jglob += nNodes_v; jcomp++;
                        }
                    } else
                    {
                        coefMat[pos] = 0;
                        symmPosToZero[countPosToZero] = pos;
                        countPosToZero++;
                    }
                }       
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            statBCAlreadyApplied = true;
        } else
        {
            //--- Wall
            int nToZero = wallPosToZero.length; int nToOne = wallPosToOne.length;
            std::shared_ptr<prec[]> coefMat = SYSMAT_NS.coef;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int iwall = 0; iwall < nWall; iwall++)
                {
                    int iglob = wallNod[iwall] + icomp*nNodes_v;
                    rhs[iglob] = 0.0;
                }
            }
            //---
            for (int countPosToZero = 0; countPosToZero < nToZero; countPosToZero++)
            {
                int pos = wallPosToZero[countPosToZero];
                coefMat[pos] = 0.0;
            }
            for (int countPosToOne = 0; countPosToOne < nToOne; countPosToOne++)
            {
                int pos = wallPosToOne[countPosToOne];
                coefMat[pos] = 1.0;
            }

            //--- NVOB
            nToZero = NVOBPosToZero.length; nToOne = NVOBPosToOne.length;
            for (int icomp = 0; icomp < dim-1; icomp++)
            {
                for (int iNVOB = 0; iNVOB < nNVOB; iNVOB++)
                {
                    int iglob = NVOBNod[iNVOB] + icomp*nNodes_v;
                    rhs[iglob] = 0.0;
                }
            }
            //---
            for (int countPosToZero = 0; countPosToZero < nToZero; countPosToZero++)
            {
                int pos = NVOBPosToZero[countPosToZero];
                coefMat[pos] = 0.0;
            }
            for (int countPosToOne = 0; countPosToOne < nToOne; countPosToOne++)
            {
                int pos = NVOBPosToOne[countPosToOne];
                coefMat[pos] = 1.0;
            }

            //--- Dirichlet
            nToOne = dirPosToOne.length; int countIglobToZero = 0;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int idir = 0; idir < nDir; idir++)
                {
                    int tempDirNod = dirNod[idir] + icomp*nNodes_v;
                    prec tempVal = dirVal[idir][icomp]; 
                    rhs[tempDirNod] = tempVal;
                    int tempDirCountId = icomp*nDir + idir;
                    int start = dirCountToZero[tempDirCountId]; int end = dirCountToZero[tempDirCountId+1];
                    int length = end - start;
                    for (int i = 0; i < length/2; i++)
                    {
                        int countPosToZero = start + i;
                        int pos = dirPosToZero[countPosToZero];
                        coefMat[pos] = 0.0;
                    }
                    for (int i = length/2; i < length; i++)
                    {
                        int countPosToZero = start + i;
                        int pos = dirPosToZero[countPosToZero];
                        int iglob = dirIglobToZero[countIglobToZero]; countIglobToZero++;
                        rhs[iglob] -= coefMat[pos] * tempVal;
                        coefMat[pos] = 0.0;
                    }
                }
            }
            //---
            for (int i = 0; i < nToOne; i++)
            {
                int pos = dirPosToOne[i];
                coefMat[pos] = 1;
            }
            //----------------------------------
            // SYMMETRY BC
            //----------------------------------
            int countIglobToCoef = 0; nToZero = symmPosToZero.length;

            for (int isymm = 0; isymm < nSymm; isymm++)
            {
                int trueIglob = symmRealIglob[isymm];

                prec* tempVal = symmVal[isymm]; 

                rhs[trueIglob] = 0;

                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    int pos = symmPosToCoef[countIglobToCoef];
                    coefMat[pos] = tempVal[jcomp];
                    countIglobToCoef++;
                }
            }
            //---
            for (int i = 0; i < nToZero; i++)
            {
                int pos = symmPosToZero[i];
                coefMat[pos] = 0;
            }
        }
        //---
        break;
    }
    }
    // NEUMANN BC
    //static
    for (int icomp = 0; icomp < dim; icomp++)
    {
        for (int ineu = 0; ineu < nNeu; ineu++)
        {
            int iglob = neuNod[ineu] + icomp*nNodes_v;
            rhs[iglob] -= neuVal[ineu][icomp];
        }
    }
}

//-----------------------------
// IMPOSE TIME DEPENDENT BC
//-----------------------------
void ADJOINT_NS::imposeTimeBC(CSRMAT &SYSMAT_NS, VECTOR &rhs)
{
    static bool timeBCAlreadyApplied = false;
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS IMPOSE TIME BC  |--\n----------\n";
    int dim = (*physics).dim;
    int flagBC = (*physics).flagBC;
    int nNodes_v = (*physics).nNodes_v;
    // DIRICHLET BC
    switch (flagBC)
    {
    case 0:
    {
        prec penalty = 1e13;
        // time dep dir
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int idir = 0; idir < nTimeDir; idir++)
            {
                int iglob = dirTimeNod[idir] + nNodes_v*icomp;  
                SYSMAT_NS(iglob,iglob,penalty);
                rhs[iglob] = penalty * dirTimeVal[idir][icomp];
            }
        }
        break;
    }
    case 1:
    {
        // time dep
        // set to zero row % col entries
        if (!timeBCAlreadyApplied)
        {
            int countPosToZero = 0; int countIglobToZero = 0; int countPosToOne = 0;
            std::shared_ptr<prec[]> coefMat = SYSMAT_NS.coef;
            dirTimeCountToZero[0] = 0;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int idir = 0; idir < nTimeDir; idir++)
                {
                    prec tempDirVal = dirTimeVal[idir][icomp];
                    int tempDirTimeNod  = dirTimeNod[idir] + icomp*nNodes_v;
                    rhs[tempDirTimeNod] = tempDirVal;
                    VECTOR_INT rowInColToZero = SYSMAT_NS.setColInRowToZero(tempDirTimeNod, dirTimePosToZero, countPosToZero, dirTimePosToOne, countPosToOne); //since the matrix has symmetric non zero entries we use the                                                                    //non 
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int iglob = rowInColToZero[i];
                        int pos = SYSMAT_NS.getPos(iglob, tempDirTimeNod);
                        prec tempVal = coefMat[pos];
                        rhs[iglob] -= tempVal * tempDirVal;
                        coefMat[pos] = 0.0;
                        dirTimeIglobToZero[countIglobToZero] = iglob; countIglobToZero++;
                        dirTimePosToZero[countPosToZero] = pos; countPosToZero++;

                        int tempDirTimeCountId = icomp*nTimeDir + idir + 1;
                        dirTimeCountToZero[tempDirTimeCountId] = countPosToZero;
                    }
                }
            }
            timeBCAlreadyApplied = true;
        } else
        {
            int nToOne = dirTimePosToOne.length;
            std::shared_ptr<prec[]> coefMat = SYSMAT_NS.coef; int countIglobToZero = 0;
            //---
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int idir = 0; idir < nTimeDir; idir++)
                {
                    int tempDirTimeNod = dirTimeNod[idir] + icomp*nNodes_v;
                    prec tempVal = dirTimeVal[idir][icomp];
                    rhs[tempDirTimeNod] = tempVal;

                    int tempDirTimeCountId = icomp*nTimeDir + idir;
                    int start = dirTimeCountToZero[tempDirTimeCountId]; int end = dirTimeCountToZero[tempDirTimeCountId+1];
                    int length = end - start;
                    for (int i = 0; i < length/2; i++)
                    {
                        int countPosToZero = start + i;
                        int pos = dirTimePosToZero[countPosToZero];
                        coefMat[pos] = 0.0;
                    }
                    for (int i = length/2; i < length; i++)
                    {
                        int countPosToZero = start + i;
                        int pos = dirTimePosToZero[countPosToZero];
                        int iglob = dirTimeIglobToZero[countIglobToZero]; countIglobToZero++;
                        rhs[iglob] -= coefMat[pos] * tempVal;
                        coefMat[pos] = 0.0;
                    }
                }
            }
            //---
            for (int i = 0; i < nToOne; i++)
            {
                int pos = dirTimePosToOne[i];
                coefMat[pos] = 1;
            }
        }
        
        //---       
        break;
    }
    }

    // NEUMANN BC
    // time dependent
    for (int icomp = 0; icomp < dim; icomp++)
    {
        for (int ineu = 0; ineu < nTimeNeu; ineu++)
        {
            int iglob = neuTimeNod[ineu] + nNodes_v*icomp;
            rhs[iglob] -= neuTimeVal[ineu][icomp];
        }
    }
}

//-------------------------------------------------
// PREPARE SOLVER
//-------------------------------------------------
void ADJOINT_NS::prepareSolver(PROBLEM_NS &NS)
{
    H.refInitialize(NS.H);
    M.refInitialize(NS.M);
    B.resize((*physics).dim);
    for (int icomp = 0; icomp < (*physics).dim; icomp++)  
    {
        B[icomp].refInitialize(NS.B[icomp]);
    }
    //
    createSysmatBase(H, B);

    nTermB = NS.nTermB;
    nEvals = NS.nEvals;
    realPos = NS.realPos;
    assembledAdgu();
    
    for (int icomp = 0; icomp < (*physics).dim; icomp++) B[icomp].dlt();
    // prepare BC update in case of lifting functions imposition
    int dim = (*physics).dim;
    if ((*physics).flagBC == 1)
    {
        std::shared_ptr<int[]> iat = SYSMAT_base.iat;
        int countPosToZero = 0;
        //--- wall
        for (int iwall = 0; iwall < nWall; iwall++)
        {
            int iglob = wallNod[iwall];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nWall;
        countPosToZero *= 2;
        wallPosToZero.initialize(countPosToZero*dim);
        wallPosToOne.initialize(nWall*dim);
        //--- static
        countPosToZero = 0;
        for (int idir = 0; idir < nDir; idir++)
        {
            int iglob = dirNod[idir];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nDir;
        dirIglobToZero.initialize(countPosToZero * dim);
        countPosToZero *= 2;
        dirPosToZero.initialize(countPosToZero * dim);
        dirPosToOne.initialize(nDir * dim);
        dirCountToZero.initialize(nDir*dim + 1);
        //--- time-dependent
        countPosToZero = 0;
        for (int idir = 0; idir < nTimeDir; idir++)
        {
            int iglob = dirTimeNod[idir];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nTimeDir;
        dirTimeIglobToZero.initialize(countPosToZero * dim);
        countPosToZero *= 2;
        dirTimePosToZero.initialize(countPosToZero * dim);
        dirTimePosToOne.initialize(nTimeDir * dim);
        dirTimeCountToZero.initialize(nTimeDir*dim + 1);

        //--- NVOB
        countPosToZero = 0;
        for (int iNVOB = 0; iNVOB < nNVOB; iNVOB++)
        {
            int iglob = NVOBNod[iNVOB];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nNVOB;
        // NVOBIglobToZero.initialize(countPosToZero * (dim-1));
        countPosToZero *= 2;
        NVOBPosToZero.initialize(countPosToZero * (dim-1));
        NVOBPosToOne.initialize(nNVOB * (dim-1));
        // NVOBCountToZero.initialize(nNVOB*(dim-1) + 1);


        //--------------------------------------------------
        // SYMMETRY BC
        //--------------------------------------------------
        countPosToZero = 0;
        for (int isymm = 0; isymm < nSymm; isymm++)
        {
            int iglob = symmNod[isymm];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nSymm*dim;
        symmPosToZero.initialize(countPosToZero);
        symmPosToCoef.initialize(nSymm*dim);
    }
    //-----------------------
    // FREE USELESS STUFF
    //-----------------------
    // dirBound.dlt();
    // wallBound.dlt(); 

    setInitCond();
}

//-------------------------------------------------
// forcing
//-------------------------------------------------
void ADJOINT_NS::applyStatForcing()
{
    int dim = (*physics).dim; int nNodes_v = (*physics).nNodes_v;
    int vecLength = dim*nNodes_v + (*physics).nNodes;

    rhs.initialize(vecLength);

}


//-------------------------------------------------
// UPDATE SYSMAT_ADJ
//-------------------------------------------------
void ADJOINT_NS::updateSYSMAT(VECTOR &lastSolNS)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| UPDATE SYSMAT_ADJ |--\n----------\n";
    SYSMAT = SYSMAT_base; // initialize SYSMAT as SYSMAT_base in order to add updated matrices in correct blocks

    // add Time Derivative Matrix
    VECTOR coefM(M.coef, nTerms);

    std::vector<std::vector<VECTOR>> coefMVec = {{coefM}};
    prec factorM = 1/deltaT;
    // add P
    addToSysmat(coefMVec, 0, factorM); // sum only on diagonal blocks
    // // //---

    // add Optimization Variable Matrix
    // assemble
    assembleMa();
    std::vector<std::vector<VECTOR>> coefMaVec = {{Ma}};
    prec factorMa = 1 / (*physics).rho;
    // add Ma
    addToSysmat(coefMaVec, 0, factorMa); // sum only on diagonal blocks
    // //---

    // add Adjoint Matrices
    // assemble
    if ((*physics).isNS == 1)
    {
        assembleA(lastSolNS); // assebmle Ac & Ag coefficients
        //add Ac
        std::vector<std::vector<VECTOR>> coefAcVec = {{Ac}};
        addToSysmat(coefAcVec); // sum only on diagonal blocks
        // add Ag
        addToSysmat(Ag, 1); // sum on all the blocks of H pattern in SYSMAT_ADJ
    }
    //---
}
//---------------------------
// SOLVER
//---------------------------
void ADJOINT_NS::StatSolver(VECTOR &lastSolNS)
{
     // pos dal primo BC
        globIter = 1;
        prec t_end = (*physics).t_end;
        deltaT = (*physics).deltaT;
        time = t_end;
        
        oneStepSolverAdjoint(lastSolNS);
        // std::cout << "\nADJ_las_sol_norm: " << lastSol.norm() << "\n";
        // pause();
        (*physics).ADJ_solution.modifyRow(0, lastSol);
        print_sol_in_VTK((*physics).ADJ_solution);
        VTKWriter.closeTFile();
}
//---
void ADJOINT_NS::Solver()
{
    int maxTimeIter = (*physics).solution_deltaT.length;
    globIter = maxTimeIter;
    // (*physics).ADJ_solution.complete_reset();
    // (*physics).ADJ_solution.initialize((*physics).NS_solution.nRow, (*physics).NS_solution.nCol);
    prec t_end = (*physics).t_end;

    time = t_end;

    lastSol.resetZeros();   // in the adj system there is no  convergence iterations,
                            // so the lastSol is not used in the solution procedure
                            // and can be reinitialized without slowing the solution procedure

    printf("|-| IT: 0 \tTIME: %7.4" format " \tFINAL CONDITION\n", t_end);
    curr_iter = maxTimeIter;

    (*physics).ADJ_solution.modifyRow(curr_iter, lastSol);
    // print_one_step_sol_in_VTK((*physics).dim, (*physics).nNodes, (*physics).nNodes_v, t_end);
    
    for (int iter = maxTimeIter-1; iter > -1; iter--)
    {
        curr_iter = iter;
        VECTOR temp_NS_sol = (*physics).NS_solution.get_row(iter);
        deltaT = (*physics).solution_deltaT[iter];
        oneStepSolverAdjoint(temp_NS_sol);
        // std::cout << "\nADJ_las_sol_norm: " << lastSol.norm() << "\n";
        (*physics).ADJ_solution.modifyRow(curr_iter, lastSol);
    }
    // pause();
    print_sol_in_VTK((*physics).ADJ_solution);
    VTKWriter.closeTFile();
}

//-------------------------
// ONE-STEP SOLVER
//-------------------------
void ADJOINT_NS::oneStepSolverAdjoint(VECTOR &nsSol)
{
    prec trialTime = time - deltaT;
    
    // updateBC(trialTime);
    updateSYSMAT(nsSol);
    updateRHS(nsSol);

    // if (completeLog) 
    if (completeLog < 2) std::cout << "\n----------------------------------------------------------------------------\n--| ADJ ONE-STEP SOLVER TIME " << time << "-> " << trialTime << "|--\n----------------------------------------------------------------------------\n";

    VECTOR rhs_final;
    CSRMAT SYSMAT_final;

    rhs_final = rhs;
    SYSMAT_final = SYSMAT;
    
    imposeBC(SYSMAT_final, rhs_final);

    for (int i = 0; i < nDof; i++) 
    {
        if (!(abs(rhs_final[i]) < 1e16)) 
        {
            std::string err = "\ni:" + std::to_string(i) + "\n";
            throw_line(err); 
        }
    }

    // call gmres
    VECTOR sol(nDof);
    std::shared_ptr<prec[]>  solP = sol.P;

    if (completeLog < 2)  std::cout << "\n----------\n--| SOLVING ADJOINT PROBLEM |--\n----------\n";

    double startTime = omp_get_wtime();
    
    //SOLVERLS::launchPardiso(SYSMAT_final, rhs_final.P, solP);

    //precondSchurSolver(SYSMAT_final, rhs_final, sol);

    prec finalRes;
    SIMPLESolver(SYSMAT_final, rhs_final, sol, finalRes);

    for (int i = 0; i < nDof; i++) 
    {
        // printf("sol[i]: %" format "\n", sol[i]);
        if (!(abs(sol[i]) < 1e16)) throw_line("eccolo 2222\n"); 
    }

    // prec bNorm   = rhs_final.norm();
    // sol = SOLVERLS::gmres_p(SYSMAT_final, rhs_final.P, lastSol.P, bNorm, 1e-6, 500); 

    // t2 = clock();
    // gettimeofday(&tv2, NULL);
    // if (completeLog < 2)  std::cout << "in " << (double)(tv2.tv_sec+tv2.tv_usec*1e-6 - (tv1.tv_sec+tv1.tv_usec*1e-6)) << " seconds (user time)\n";
    // if (completeLog < 2)  std::cout << "in " << (double)(t2-t1) / CLOCKS_PER_SEC << " seconds (CPU time)\n";

    double endTime = omp_get_wtime();

    switch ((*physics).isNS)
    {
        case 0:
        {
            printf("|ADJ_S| IT: %d \tTIME: %7.4" format " \tSOLVER TIME: %" format "\n", ((*physics).solution_deltaT.length - globIter + 1), trialTime, endTime-startTime);
            break;
        }
        case 1:
        {
            printf("|ADJ_NS| IT: %d \tTIME: %7.4" format " \tSOLVER TIME: %" format "\n", ((*physics).solution_deltaT.length - globIter + 1), trialTime, endTime-startTime);
            break;
        }
        default:
        {
            throw_line("ERROR: not handled problem physics (Fluid physics different from Stokes or Navier-Stokes)\n");
            break;
        }
    }
    // printf("|ADJ| IT: %d \tTIME: %7.4" format " \tSOLVER TIME: %" format "\n", ((*physics).solution_deltaT.length - globIter + 1), trialTime, endTime-startTime);
    
    PARALLEL::copy(sol, lastSol);

    time = trialTime;
    //solution_times.append(time);
    // SYSMAT_final.dlt();

    // // GET VELOCITY
    // MATRIX velocity = getVelocityFromSol2(lastSol);
    // // GET PRESSURE
    // VECTOR pressure(nNodes);
    // prec* tempPointer = &lastSol[dim*nNodes_v];
    // pressure = tempPointer;
    // std::vector<VECTOR> analyticSolV(dim);
    // //-----------------
    // if (printRes) VTKWriter.write((*physics).coord, (*physics).elem, time, pressure, 1, "Pressure", velocity, dim, "Velocity");
    
    // (*physics).ADJ_solution.modifyRow(curr_iter, lastSol);

    // print_one_step_sol_in_VTK(dim, nNodes, nNodes_v, time);
    
    globIter--;
    
    // std::ofstream resFile;
    // std::string resFilePath = "ADJCurrSol/" + std::to_string(globIter) + ".txt";
    // resFile.open(resFilePath, std::ios::out | std::ios::binary);
    // if(!resFile) {
    //   throw_line("ERROR: Can't open NS result file\n");
    // }
    // for (int i = 0; i < nDof; i++) resFile.write((char *) &lastSol[i], sizeof(prec));
    // resFile.close();
}
//----------------------------------------------------


//-------------------------------------------------
// UPDATE RHS
//-------------------------------------------------
void ADJOINT_NS::updateRHS(VECTOR &nsSol)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| UPDATE RHS |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v; int nDof_v = dim*nNodes_v;

    // updateTimeForcing();
    int nElem_v = (*physics).nElem_v;
    VECTOR Volume_v = (*physics).Volume_v;
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);

    prec rho = (*physics).rho;
    // prec mu  = (*physics).mu;
    int vecLength = dim*nNodes_v + (*physics).nNodes;

    rhs.initialize(vecLength);
    PARALLEL::resetZeros(rhs);

    VECTOR lastSolU(nDof_v);
    prec* nsSolP = &(nsSol[0]); 
    prec* rhsP = &(rhs[0]);
    PARALLEL::prod(dAdGu, nsSolP, nDof_v, rhsP); //rhs *= -1;

    for (int i = 0; i < nDof; i++) 
    {
        if (!(abs(rhs[i]) < 1e16)) 
        {
            std::string err = "\ni:" + std::to_string(i) + "\n";
            throw_line(err); 
        }
    }

    //--------------------
    for (int icomp = 0; icomp < dim; icomp++)
    {
        VECTOR UcompSol(nNodes_v);
        for (int i = 0; i < nNodes_v; i++) UcompSol[i] = nsSol[icomp*nNodes_v+i];

        for (int iel = 0; iel < nElem_v; iel++)
        {
            int* tempElem = elem_v[iel];

            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = tempElem[iloc];

                for (int jloc = 0; jloc < dim+1; jloc++)
                {
                    int countEqual = 0;
                    int jglob = tempElem[jloc];

                    prec tempCompVel = UcompSol[jglob];

                    if (iloc == jloc) countEqual++;

                    for (int kloc = 0; kloc < dim+1; kloc++)
                    {
                        int kglob = tempElem[kloc];
                        prec factor = 1;

                        int tempCount = countEqual;

                        if (kloc == iloc) tempCount++; 
                        if (kloc == jloc) tempCount++;

                        switch (tempCount)
                        {
                            case 0:
                            {
                                factor = 1.0/120;
                                break;
                            }
                            case 1:
                            {
                                factor = 1.0/60;
                                break;
                            }
                            case 3: 
                            {
                                factor = 1.0/20;
                            }
                        }

                        factor *= (4-dim)*Volume_v[iel];
                        rhs[iglob+nNodes_v*icomp] -= (2*fWeights[0]*factor*alpha[kglob]/rho) *tempCompVel; // fWeights[4] * tempCompVel
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < nDof; i++) 
    {
        if (!(abs(rhs[i]) < 1e16)) 
        {
            std::string err = "\ni:" + std::to_string(i) + "\n";
            throw_line(err); 
        }
    }
    
    rhs /= (*physics).func_normalization_factor;

    // functional constraints contribution
    for (int icons = 0; icons < constraints->n_constr; icons++)
    {    
        int type = constraints->list[icons].type;
        switch (type)
        {
            case 4:
            {
                update_WSS_constraint(g, icons, constraints.list[icons]);
                break;
            }
            default:
            {
                break;
            }  
        }
    }

    // time iteration contribution
    for (int icomp = 0; icomp < dim; icomp++)
    {
        int startId = icomp*nNodes_v;
        VECTOR tempSolComp(nNodes_v);
        prec* tempP = &lastSol[startId];
        tempSolComp.copy(tempP);

        VECTOR tSol(nNodes_v);
        M.prod(tempSolComp, tSol);
        // tempSolComp = M*tempSolComp;
        tSol /= deltaT;
        tempP = &rhs[startId];
        VECTOR::sum(tempP, tSol.P, nNodes_v, tempP);
    }

    for (int i = 0; i < nDof; i++) 
    {
        if (!(abs(rhs[i]) < 1e16)) 
        {
            std::string err = "\ni:" + std::to_string(i) + "\n";
            throw_line(err); 
        }
    }
    
}