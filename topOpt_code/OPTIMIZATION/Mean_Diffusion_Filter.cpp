#include "Mean_Diffusion_Filter.h"

void Mean_DIFFUSION_FILTER::initialize(int diff_filter_case, PHYSICS* physicsP, int nNodesInDomain, VECTOR_INT nodesInDomain, VECTOR_INT optimizationNodeFromGlobNode, MATRIX topologyOptBox, prec diffRadiusPercentace, prec diffusionFilterWeight)
    {
        physics = physicsP;
        dim = (*physics).dim;
        nNodesInDom = nNodesInDomain;
        nodesInDom = nodesInDomain;
        diffusion_filter_case = diff_filter_case;
        optNodeFromGlobNode = optimizationNodeFromGlobNode;
        topOptBox = topologyOptBox;
        topOptBoxSize.initialize(dim);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            topOptBox[icomp][0] -= 1e-10;
            topOptBox[icomp][1] += 1e-10;
            topOptBoxSize[icomp] = abs(topOptBox[icomp][1] - topOptBox[icomp][0]);
        }
        diffRadius = diffRadiusPercentace * topOptBoxSize.min();
        cellSize = 2*diffRadius;
        diffFilterWeight = diffusionFilterWeight;
        initializeDiffFilter();
    }

void Mean_DIFFUSION_FILTER::initializeDiffFilter()
{
    nCells.setZeros(3);
    nCells += 1;
    for (int icomp = 0; icomp < dim; icomp++)
    {
        nCells[icomp] = ceil(topOptBoxSize[icomp] / cellSize);
    }
    cellNodesTensor.resize(nCells[0]);
    for (int icell = 0; icell < nCells[0]; icell++)
    {
        cellNodesTensor[icell].resize(nCells[1]);
        for (int jcell = 0; jcell < nCells[1]; jcell++)
        {
            cellNodesTensor[icell][jcell].resize(nCells[2]);
            for (int kcell = 0; kcell < nCells[2]; kcell++)
            {
                cellNodesTensor[icell][jcell][kcell].initialize(0);
            }                  
        }
    }

    buildCellsTensor();
    buildNodesNB();
}

void Mean_DIFFUSION_FILTER::buildCellsTensor()
{
    MATRIX* coord_v = &((*physics).coord_v);
    for (int inode = 0; inode < nNodesInDom; inode++)
    {
        int iglob = nodesInDom[inode];
        prec* tempCoord = (*coord_v)[iglob];
        VECTOR_INT cellIds;
        cellIds.setZeros(3);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            cellIds[icomp] = floor((tempCoord[icomp] - topOptBox[icomp][0]) / cellSize);
        }
        cellNodesTensor[cellIds[0]][cellIds[1]][cellIds[2]].append(iglob);
    }
}

void Mean_DIFFUSION_FILTER::buildNodesNB()
{
    nodesNB.resize(nNodesInDom);
    VECTOR_INT localNB(3);
    localNB.set3Values(-1, 0, 1);

    for (int icell = 0; icell < nCells[0]; icell++)
    {
        for (int jcell = 0; jcell < nCells[1]; jcell++)
        {
            for (int kcell = 0; kcell < nCells[2]; kcell++)
            {
                VECTOR_INT possibleNB(0);
                for (int iloc = 0; iloc < localNB.length; iloc ++)
                {
                    int icellNB = icell + localNB[iloc];
                    if ((icellNB > -1 ) && (icellNB < nCells[0]))
                    {
                        for (int jloc = 0; jloc < localNB.length; jloc ++)
                        {
                            int jcellNB = jcell + localNB[jloc];
                            if ((jcellNB > -1 ) && (jcellNB < nCells[1]))
                            {
                                for (int kloc = 0; kloc < localNB.length; kloc ++)
                                {
                                    int kcellNB = kcell + localNB[kloc];
                                    if ((kcellNB > -1 ) && (kcellNB < nCells[2]))
                                    {
                                        VECTOR_INT tempNB = cellNodesTensor[icellNB][jcellNB][kcellNB];
                                        possibleNB.append(tempNB);
                                    }
                                }  
                            }
                        }
                    }
                }

                VECTOR_INT tempCellNodes = cellNodesTensor[icell][jcell][kcell];
                for (int inode = 0; inode < tempCellNodes.length; inode++)
                {
                    int tempGlobNode = tempCellNodes[inode];
                    int tempOptNode = optNodeFromGlobNode[tempGlobNode];
                    nodesNB[tempOptNode].initialize(dim, diffRadius, diffFilterWeight, tempGlobNode, tempOptNode, (*physics).coord_v[tempGlobNode], possibleNB);
                }
            }
        }
    }

    buildNeighbourhoods();
}

void Mean_DIFFUSION_FILTER::buildNeighbourhoods()
{
    std::cout << "\nBUILD NEIGHBOURHOODS \n";
    
    int node_counter = 0;
    #pragma omp parallel num_threads(std::thread::hardware_concurrency())
    {
        #pragma omp for
        for (int inode = 0; inode < nNodesInDom; inode++)
        {
            #pragma omp critical (Build_NB_progress_bar)
            {
                node_counter++;
                std::cout << "--| Build NB | Node: " << node_counter << "/" << nNodesInDom << "\t\tperc:" << floor((double(node_counter) / double(nNodesInDom))*100) << "\n";
            }
            nodesNB[inode].buildNeighbourhood_v1(nodesNB, optNodeFromGlobNode);
        }
    }
}
    

void Mean_DIFFUSION_FILTER::filterGamma_v0(VECTOR &gammaOld, VECTOR &gammaNew)
{
    gammaNew = gammaOld;
    for (int inode = 0; inode < nNodesInDom; inode++)
    {
        Mean_DF_NODE_NB nodeNB = nodesNB[inode];
        int iglob = nodeNB.globId;
        if (nodeNB.neighbours.length != 0)
        {
            prec gammaLocNBAvg = VECTOR::weigthedMeanByIdSet(gammaOld, nodeNB.neighbours, nodeNB.weigthNB);
            gammaNew[iglob] = (1-diffFilterWeight) * gammaOld[iglob] + diffFilterWeight* gammaLocNBAvg;
        }
    }
}

void Mean_DIFFUSION_FILTER::filter_gamma_and_derivative(VECTOR &gamma, VECTOR &gamma_filter, VECTOR &dgamma_filter)
{
    filter_gamma(gamma, gamma_filter, diffusion_filter_case);
    eval_gamma_filter_derivative(gamma, gamma_filter, dgamma_filter);
}

void Mean_DIFFUSION_FILTER::filter_gamma(VECTOR &gamma, VECTOR &gamma_filter, int filter_case)
{
    for (int inode = 0; inode < nNodesInDom; inode++)
    {
        Mean_DF_NODE_NB* nodeNB = &(nodesNB[inode]);
        //int iglob = (*nodeNB).globId;
        if ((*nodeNB).neighbours.length == 0) 
        {
            throw_line("ERROR in gamma filtering: the node itself isn't in its neighbourood\n");
        }
        else
        {
            switch (filter_case)
            {
                case 1:
                {
                    prec gammaLocNBAvg = VECTOR::weigthedMeanByIdSet(gamma, (*nodeNB).neighbours, (*nodeNB).weigthNB);
                    gamma_filter[inode] = gammaLocNBAvg;
                    break;
                }
                case 2:
                {
                    prec gammaLocNBAvg = VECTOR::geometricMeanByIdSet(gamma, (*nodeNB).neighbours, (*nodeNB).weigthNB);
                    gamma_filter[inode] = gammaLocNBAvg;
                    break;
                }                
                default:
                {
                    throw_line("ERROR: evaluating a Mean based filter on gamma with a wrong diffusive filter case value.\n");
                    break;
                }
            }
        }
    }
}

void Mean_DIFFUSION_FILTER::eval_gamma_filter_derivative(VECTOR &gamma, VECTOR &gamma_filter, VECTOR &dgamma_filter)
{
    for (int inode = 0; inode < nNodesInDom; inode++)
            {
                dgamma_filter[inode] = nodesNB[inode].eval_gamma_derivative_weight(diffusion_filter_case, gamma, gamma_filter);
            }
}

void Mean_DIFFUSION_FILTER::printNeighbourhood()
{
    for (int inode = 0; inode < nNodesInDom; inode++)
    {
        nodesNB[inode].printNeighbours();
        std::cout << "\n\n";
    }
}