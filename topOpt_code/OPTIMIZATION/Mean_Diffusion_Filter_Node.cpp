#include "Mean_Diffusion_Filter.h"

void Mean_DF_NODE_NB::initialize(int dimension, prec diffusionRadius, prec diffusionFilterWeight, int nodeGlobId, int nodeOptId, prec* coords, VECTOR_INT possibleNeighbours)
{
    dim = dimension;
    diffRadius = diffusionRadius;
    diffFilterWeight = diffusionFilterWeight;
    globId = nodeGlobId;
    optId = nodeOptId;
    coord.initialize(coords, dim);
    possibleNB = possibleNeighbours;
}

void Mean_DF_NODE_NB::buildNeighbourhood_v0(std::vector<Mean_DF_NODE_NB> &nodesNB, VECTOR_INT &optNodeFromGlobNode)
{
    int nMaxNB = possibleNB.length;
    neighbours.initializeZero(nMaxNB);
    VECTOR distNB(nMaxNB);
    prec weigthsSum = 0.0;
    weigthNB.setZeros(nMaxNB);
    int nNB = 0;
    for (int inb = 0; inb < nMaxNB; inb++)
    {
        int idNB = possibleNB[inb];
        int optIdNB = optNodeFromGlobNode[idNB];
        Mean_DF_NODE_NB tempNode = nodesNB[optIdNB];
        if (globId != tempNode.globId)
        {
            VECTOR coordNB = tempNode.coord;
            //std::cout << "Pre diff\n";
            if (coordNB.length != dim) 
            {
                coord.printRowMatlab("Coord");
                coordNB.printRowMatlab("CoordNB");
                std::string errorStr = "ERROR inserting a neighbour for the opt node " + std::to_string(inb) + ". Neighbour optId:" + std::to_string(optIdNB) + "\n";
                throw_line(errorStr);
            }
            VECTOR distVec = coord - coordNB;
            //std::cout << "Post diff\n";
            prec tempDist = VECTOR::norm(distVec);
            //std::cout << "Dist: " << tempDist << "; Radius: " << diffRadius << "; CHECK: " << (tempDist <= diffRadius) << "\n";
            if (tempDist <= diffRadius)
            {
                neighbours[nNB] = idNB;
                prec tempWeigth = 1 / tempDist;
                weigthNB[nNB] = tempWeigth;
                weigthsSum += tempWeigth;
                nNB += 1;
            }
        }
    }
    neighbours.shrink(nNB);
    weigthNB.shrink(nNB);
    weigthNB = weigthNB / weigthsSum;
}

void Mean_DF_NODE_NB::buildNeighbourhood_v1(std::vector<Mean_DF_NODE_NB> &nodesNB, VECTOR_INT &optNodeFromGlobNode)
{
    int nMaxNB = possibleNB.length;
    neighbours.initializeZero(nMaxNB);
    VECTOR distNB(nMaxNB);
    prec weigthsSum = 0.0;
    weigthNB.setZeros(nMaxNB);
    int nNB = 0;
    for (int inb = 0; inb < nMaxNB; inb++)
    {
        int idNB = possibleNB[inb];
        int optIdNB = optNodeFromGlobNode[idNB];
        Mean_DF_NODE_NB tempNB = nodesNB[optIdNB];
        VECTOR coordNB = tempNB.coord;
        if (coordNB.length != dim) 
        {
            #pragma omp critical (error_in_coordNB)
            {
                coord.printRowMatlab("Coord");
                coordNB.printRowMatlab("CoordNB");
                std::string errorStr = "ERROR inserting a neighbour for the opt node " + std::to_string(inb) + ". Neighbour optId:" + std::to_string(optIdNB) + "\n";
                throw_line(errorStr);
            }
            
        }
        VECTOR distVec = coord - coordNB;
        prec tempDist = VECTOR::norm(distVec);
        if (tempDist <= diffRadius)
        {
            neighbours[nNB] = optIdNB;
            prec tempWeigth = tempDist / diffRadius;
            tempWeigth = 1 - tempWeigth;
            if (globId != tempNB.globId)
            {
                tempWeigth *= diffFilterWeight;
            }
            weigthNB[nNB] = tempWeigth;
            weigthsSum += tempWeigth;
            nNB += 1;
        }
    }
    neighbours.shrink(nNB);
    weigthNB.shrink(nNB);
    weigthNB = weigthNB / weigthsSum;

    #pragma omp critical (build_weight_as_NB)
    {
        build_weight_as_NB(nodesNB);
    }
}

void Mean_DF_NODE_NB::build_weight_as_NB(std::vector<Mean_DF_NODE_NB> &nodesNB)
{
    for (int inb = 0; inb < neighbours.length; inb++)
    {
        int optIdNB = neighbours[inb];
        Mean_DF_NODE_NB* tempNB = &(nodesNB[optIdNB]);
        (*tempNB).weight_as_NB[optId] = weigthNB[inb];
    }
}

prec Mean_DF_NODE_NB::eval_gamma_derivative_weight(int filter_case, VECTOR &gamma, VECTOR &gamma_filter)
{
    prec temp_derivative_weight = 0.0;
    prec gamma_node = gamma[optId];
    for (auto weight_pair : weight_as_NB)
    {
        prec temp_weigth_as_nb = weight_pair.second;
        if (filter_case == 2)
        {
            int node_nb = weight_pair.first;
            temp_weigth_as_nb *= gamma_filter[node_nb] / gamma_node;
        }
        temp_derivative_weight += temp_weigth_as_nb;
    }
    gamma_derivative_weight = temp_derivative_weight;
    return temp_derivative_weight;
}

void Mean_DF_NODE_NB::printNeighbours()
{
    std::string vecName = "Node" + std::to_string(globId) + " NB";
    neighbours.printRowMatlab(vecName);
}