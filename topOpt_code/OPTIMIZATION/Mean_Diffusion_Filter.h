#pragma once

#include "../CODE_HEADERS/codeHeader.h"
#include "../geometry.h"

class Mean_DF_NODE_NB //NB stands for Neighbourhood
{
public:
    int dim;
    prec diffRadius;
    prec diffFilterWeight;
    int globId;
    int optId;
    VECTOR coord;
    VECTOR_INT possibleNB;
    VECTOR_INT neighbours;
    VECTOR weigthNB;
    std::map<int, prec> weight_as_NB; //works with optmization nodes indices, contains the weight of the current node as nb of its own neighbours
    prec gamma_derivative_weight = 0.0;

    Mean_DF_NODE_NB(){};

    Mean_DF_NODE_NB(int dimension, prec diffusionRadius, prec diffusionFilterWeight, int nodeGlobId, int nodeOptId, prec* coords, VECTOR_INT possibleNeighbours)
    {
        dim = dimension;
        diffRadius = diffusionRadius;
        diffFilterWeight = diffusionFilterWeight;
        globId = nodeGlobId;
        optId = nodeOptId;
        prec* coordsP = coords;
        coord.initialize(coordsP, dim);
        possibleNB = possibleNeighbours;
    }


    void initialize(int dimension, prec diffusionRadius, prec diffusionFilterWeight, int nodeGlobId, int nodeOptId, prec* coords, VECTOR_INT possibleNeighbours);

    void buildNeighbourhood_v0(std::vector<Mean_DF_NODE_NB> &nodesNB, VECTOR_INT &optNodeFromGlobNode);
    void buildNeighbourhood_v1(std::vector<Mean_DF_NODE_NB> &nodesNB, VECTOR_INT &optNodeFromGlobNode);

    void build_weight_as_NB(std::vector<Mean_DF_NODE_NB> &nodesNB);

    prec eval_gamma_derivative_weight(int filter_case, VECTOR &gamma, VECTOR &gamma_filter);

    void printNeighbours();
};

class Mean_DIFFUSION_FILTER
{
public: 
    PHYSICS* physics;
    int dim;
    int nNodesInDom;
    int diffusion_filter_case;
    VECTOR_INT nodesInDom;
    VECTOR_INT optNodeFromGlobNode;
    MATRIX topOptBox;
    VECTOR topOptBoxSize;
    prec diffRadius;
    prec diffFilterWeight;
    prec cellSize;
    VECTOR_INT nCells;
    std::vector<std::vector<std::vector<VECTOR_INT>>> cellNodesTensor;
    std::vector<Mean_DF_NODE_NB> nodesNB;

    Mean_DIFFUSION_FILTER(){};

    Mean_DIFFUSION_FILTER(int diff_filter_case, PHYSICS* physicsP, int nNodesInDomain, VECTOR_INT nodesInDomain, VECTOR_INT optimizationNodeFromGlobNode, MATRIX topologyOptBox, prec diffRadiusPercentace, prec diffusionFilterWeight)
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
            topOptBoxSize[icomp] = abs(topOptBox[icomp][1] - topOptBox[icomp][0]);
        }
        diffRadius = diffRadiusPercentace * topOptBoxSize.min();
        cellSize = 2*diffRadius;
        diffFilterWeight = diffusionFilterWeight;
        initializeDiffFilter();
    }

    void initialize(){};

    void initialize(int diff_filter_case, PHYSICS* physicsP, int nNodesInDomain, VECTOR_INT nodesInDomain, VECTOR_INT optimizationNodeFromGlobNode, MATRIX topologyOptBox, prec diffRadiusPercentace, prec diffusionFilterWeight);

    void initializeDiffFilter();

    void buildCellsTensor();

    void buildNodesNB();

    void buildNeighbourhoods();

    void filter_gamma_and_derivative(VECTOR &gamma, VECTOR &gamma_filter, VECTOR &dgamma_filter);

    void filterGamma_v0(VECTOR &gamma, VECTOR &gamma_filter);

    void filter_gamma(VECTOR &gamma, VECTOR &gamma_filter, int filter_case);

    void eval_gamma_filter_derivative(VECTOR &gamma, VECTOR &gamma_filter, VECTOR &dgamma_filter);

    void printNeighbourhood();
};
//---