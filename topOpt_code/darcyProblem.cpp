#include "CODE_HEADERS/codeHeader.h"
#include "problemDarcyHeader.h"


//----------------------------------------------------
// INITIALIZE a DARCY PROBLEM
//----------------------------------------------------
void PROBLEM_DARCY::initialize(PHYSICS *&Physics, std::string probRefFile, VECTOR &alphaIn, bool print)
    {
        printRes = print;
        physics = Physics;
        importParameters(probRefFile);
        checkImportParameters();
        importPREPRO();
        if (abs(time) < 1e-16) time = 0;
        localBasis();
    }

//----------------------------------------------------
// DARCY PROBLEM import parmeters
//----------------------------------------------------
void PROBLEM_DARCY::importParameters(std::string readFile)
{
    std::cout << "\n----------\n--| DARCY IMPORT PARAMETERS |--\n----------\n";
    VECTOR_INT isRepeatedID(0);
    // BEGIN STREAMING
    std::ifstream ParameterFile;
    ParameterFile.open(readFile);
    if (!ParameterFile.is_open()) throw_line("ERROR, can't open input data file");
    std::string line;
    std::istringstream iss;

    STREAM::getLines(ParameterFile, line, 1);

    // problem name 
    STREAM::getValue(ParameterFile, line, iss, name);
    (*physics).name = name;

    STREAM::getLines(ParameterFile, line, 1); 
    
    // dimension   
    STREAM::getValue(ParameterFile, line, iss, (*physics).dim);
    int dim = (*physics).dim;
    if (dim != 2 && dim != 3) throw_line("ERROR: Inizializing a Problem of dimension different from 2 or 3 \n");

    STREAM::getLines(ParameterFile, line, 3); 

    // get time parameters
    bool isStat;
    STREAM::getValue(ParameterFile, line, iss, isStat);
    (*physics).isStationary = isStat;
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, (*physics).convergence_scale_factor);
    STREAM::getLines(ParameterFile, line, 2);  
    STREAM::getValue(ParameterFile, line, iss, (*physics).t_end);
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, deltaT);
    deltaT0 = deltaT;
    (*physics).deltaT = deltaT;
    if ((*physics).deltaT > (*physics).t_end)
    {
        throw_line("ERROR: inserted a deltaT greater than the final time\n");
    }
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, (*physics).deltaT_min);

    STREAM::getLines(ParameterFile, line, 2);
    
    // FLUID PROPERTIES: rho, mu
    STREAM::getValue(ParameterFile, line, iss, (*physics).rho);
    STREAM::getValue(ParameterFile, line, iss, (*physics).mu);
    (*physics).ni = (*physics).mu/(*physics).rho;

    // MATERIAL PROPERTIES: permeability
    STREAM::getLines(ParameterFile, line, 3);
    STREAM::getValue(ParameterFile, line, iss, n_domains);
    STREAM::getLines(ParameterFile, line, 1);
    domains_permeability_priority.initialize(n_domains);
    STREAM::getRowVector(ParameterFile, line, iss, domains_permeability_priority);
    STREAM::getLines(ParameterFile, line, 1);
    (*physics).domains_permeability.initialize(n_domains);
    STREAM::getRowVector(ParameterFile, line, iss, (*physics).domains_permeability);
    domains_permeability = (*physics).domains_permeability;

    STREAM::getLines(ParameterFile, line, 3);
    // flag Forcing
    STREAM::getValue(ParameterFile, line, iss, flagForcing);
    if (flagForcing != 0 && flagForcing != 1) throw_line("ERROR: flagForcing different from 0 or 1\n");
    // FORCING
    STREAM::getLines(ParameterFile, line, 1);
    statForcing.resize(1);
    STREAM::getValue(ParameterFile, line, iss, statForcing[0]);
    STREAM::getLines(ParameterFile, line, 1);
    timeForcing.resize(1);
    STREAM::getValue(ParameterFile, line, iss, timeForcing[0]);
    
    STREAM::getLines(ParameterFile, line, 3);

    // flag BC
    STREAM::getValue(ParameterFile, line, iss, (*physics).flagBC); 
    int flagBC = (*physics).flagBC;
    if (flagBC != 0 && flagBC != 1) throw_line("ERROR: flagBC different from 0 or 1\n");

    STREAM::getLines(ParameterFile, line, 4);
    
    // Inner Walls
    STREAM::getValue(ParameterFile, line, iss, nInnerBound);
    getline(ParameterFile, line);
    innerBound.initialize(nInnerBound);
    if (nInnerBound == 0) getline(ParameterFile, line);
    else 
    {
        int *tempInnerP = &(innerBound[0]);
        STREAM::getRowVector(ParameterFile, line, iss, tempInnerP, nInnerBound);
    }
    //-----------------------------------------
    // SYMMETRIES
    STREAM::getLines(ParameterFile, line, 4);
    STREAM::getValue(ParameterFile, line, iss, nSymmBound);
    getline(ParameterFile, line);
    symmBound.initialize(nSymmBound);
    if (nSymmBound == 0) getline(ParameterFile, line);
    else 
    {
        int *tempInnerP = &(symmBound[0]);
        STREAM::getRowVector(ParameterFile, line, iss, tempInnerP, nSymmBound);
    }

    STREAM::getLines(ParameterFile, line, 2);

    // flag noFlux    
    STREAM::getValue(ParameterFile, line, iss, flagNoFlux);
    if (flagNoFlux != 0 && flagNoFlux != 1) throw_line("ERROR: flagNoFlux different from 0 or 1\n"); 

    STREAM::getLines(ParameterFile, line, 4);

    // BOUNDARY CONDITIONS
    // # No-Flux boundaries
    STREAM::getValue(ParameterFile, line, iss, nNoFluxBound);

    getline(ParameterFile, line);

    // noFlux boundaries
    noFluxBound.initialize(nNoFluxBound);
    int* tempP = &(noFluxBound.P[0]);
    STREAM::getRowVector(ParameterFile, line, iss, tempP, nNoFluxBound);

    for (int ibound = 0; ibound < nNoFluxBound; ibound ++)
    {
        if (innerBound.hasIn(noFluxBound[ibound])) throw_line("ERROR: Repeated Inenr Bound ID in noFlux bound\n");
        if (symmBound.hasIn(noFluxBound[ibound])) throw_line("ERROR: Repeated Symmetry Bound ID in noFlux bound\n");
    }

    STREAM::getLines(ParameterFile, line, 3);
    
    // # Dirichlet boundaries
    STREAM::getValue(ParameterFile, line, iss, nDirBound);

    getline(ParameterFile, line);

    // Dirichlet boundaries
    dirBound.initialize(nDirBound);
    if (nDirBound == 0) getline(ParameterFile, line);
    else 
    {
        tempP = &(dirBound[0]);
        STREAM::getRowVector(ParameterFile, line, iss, tempP, nDirBound);
    }
    if (nDirBound != 0) getline(ParameterFile, line);
    for (int ibound = 0; ibound < nDirBound; ibound ++)
    {
        if (noFluxBound.hasIn(dirBound[ibound])) throw_line("ERROR: Repeated noFlux Bound ID in Static Dirichlet\n");
        if (innerBound.hasIn(dirBound[ibound])) throw_line("ERROR: Repeated Inner Bound ID in Static Dirichlet\n");
        if (symmBound.hasIn(dirBound[ibound])) throw_line("ERROR: Repeated Symmetry Bound ID in Static Dirichlet\n");
    }
    (*physics).inlet_bounds = dirBound;

    //----------------------------------------------
    dirFunc.resize(nDirBound);
    for (int ifunc = 0; ifunc < nDirBound; ifunc++) 
    {
        dirFunc[ifunc].resize(1);
        STREAM::getRowVector(ParameterFile, line, iss, dirFunc[ifunc]);
    }

    STREAM::getLines(ParameterFile, line, 4);

    // # Neumann Boundaries
    STREAM::getValue(ParameterFile, line, iss, nNeuBound);

    getline(ParameterFile, line);

    // Neumann boundaries
    neuBound.initialize(nNeuBound);
    if (nNeuBound != 0) 
    {
        tempP = &(neuBound[0]);
        STREAM::getRowVector(ParameterFile, line, iss, tempP, nNeuBound);
    } else getline(ParameterFile, line);
    for (int ibound = 0; ibound < nNeuBound; ibound ++)
    {
        if (noFluxBound.hasIn(neuBound[ibound])) throw_line("ERROR: Repeated noFlux Bound ID in Static Neumann\n");
        if (dirBound.hasIn(neuBound[ibound])) throw_line("ERROR: Repeated Static Dirichlet Bound ID in Static Neumann\n");
        if (innerBound.hasIn(neuBound[ibound])) throw_line("ERROR: Repeated Inner Bound ID in Static Neumann\n");
        if (symmBound.hasIn(neuBound[ibound])) throw_line("ERROR: Repeated Symmetry Bound ID in Static Neumann\n");
    }
    // Neumann flux
    neuMeanFlux.initialize(nNeuBound);
    if (nNeuBound != 0) 
    {
        getline(ParameterFile, line);
        STREAM::getColVector(ParameterFile, line, iss, neuMeanFlux, nNeuBound);
    }

    // //-------------------------------------------------
    // //------------------------------------
    // // temporal BCs
    // //------------------------------------
    STREAM::getLines(ParameterFile, line, 4);

    // Dirichlet Time Dependent
    STREAM::getValue(ParameterFile, line, iss, nDirTimeBound);
    dirTimeBound.initialize(nDirTimeBound);
    nIdDirTimeCases.initialize(nDirTimeBound);
    dirTimeFunc.resize(nDirTimeBound);
    getline(ParameterFile, line);
    for (int ibound = 0; ibound < nDirTimeBound; ibound++)
    {
        STREAM::getDirTime(ParameterFile, line, iss, 1, dirTimeBound[ibound], nIdDirTimeCases[ibound], dirTimeFunc[ibound]);
        getline(ParameterFile, line);
    }

    for (int ibound = 0; ibound < nDirTimeBound; ibound ++)
    {
        if (noFluxBound.hasIn(dirTimeBound[ibound])) throw_line("ERROR: Repeated noFlux Bound ID in Time Dependent Dirichlet\n");
        if (dirBound.hasIn(dirTimeBound[ibound])) throw_line("ERROR: Repeated Static Dirichlet Bound ID in Time Dependent Dirichlet\n");
        if (neuBound.hasIn(dirTimeBound[ibound])) throw_line("ERROR: Repeated Static Neumann Bound ID in Time Dependent Dirichlet\n");
        if (innerBound.hasIn(dirTimeBound[ibound])) throw_line("ERROR: Repeated Inner Bound ID in Time Dependent Dirichlet\n");
        if (symmBound.hasIn(dirTimeBound[ibound])) throw_line("ERROR: Repeated Symmetry Bound ID in Time Dependent Dirichlet\n");
    }
    
    STREAM::getLines(ParameterFile, line, 3);

    // Neumann Time Dependent
    STREAM::getValue(ParameterFile, line, iss, nNeuTimeBound);
    neuTimeBound.initialize(nNeuTimeBound);
    nIdNeuTimeCases.initialize(nNeuTimeBound);
    neuTimeMeanFlux.resize(nNeuTimeBound);
    getline(ParameterFile, line);
    for (int ibound = 0; ibound < nNeuTimeBound; ibound++)
    {
        STREAM::getNeuTime(ParameterFile, line, iss, neuTimeBound[ibound], nIdNeuTimeCases[ibound], neuTimeMeanFlux[ibound]);
        getline(ParameterFile, line);
    }

    for (int ibound = 0; ibound < nNeuTimeBound; ibound ++)
    {
        if (noFluxBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated noFlux Bound ID in Time Dependent Neumann\n");
        if (dirBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated Static Dirichlet Bound ID in Time Dependent Neumann\n");
        if (neuBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated Static Neumann Bound ID in Time Dependent Neumann\n");
        if (dirTimeBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated Time Dependent Dirichlet Bound ID in Time Dependent Neumann\n");
        if (innerBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated Inner Bound ID in Time Dependent Neumann\n");
        if (symmBound.hasIn(neuTimeBound[ibound])) throw_line("ERROR: Repeated Symmetry Bound ID in Time Dependent Neumann\n");
    }
    
    // CLOSE STREAMING
    ParameterFile.close();
}

//--------------------------------------------------
// IMPORT PREPRO
//--------------------------------------------------
void PROBLEM_DARCY::importPREPRO()
{
    // READ:
    // only velocity NODES
    // velocity and pressure ELEMENTS
    // only velocity BOUND INFO

    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| DARCY IMPORT PREPRO |--\n----------\n";
    std::string folderPath    = name;
    folderPath = "PREPRO/PROBLEM_DATA/" + folderPath;
    std::string NodeFilePath  = folderPath + "/Nodes.txt";
    std::string ElemFilePath  = folderPath + "/Elems.txt";
    std::string BoundNodeFilePath = folderPath + "/BoundNodes.txt";
    std::string BoundElemFilePath = folderPath + "/BoundElems.txt";
    std::string ElemGeoIdsFilePath = folderPath + "/ElemGeoIds.txt";
    (*physics).bound_info_file_path = BoundElemFilePath;

    //-------OPEN STREAMS-------
    std::ifstream  NodeFile;  NodeFile.open(NodeFilePath, std::ios::in | std::ios::binary);  if (!NodeFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  ElemFile;  ElemFile.open(ElemFilePath, std::ios::in | std::ios::binary);  if (!ElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  BoundNodeFile; BoundNodeFile.open(BoundNodeFilePath, std::ios::in | std::ios::binary); if (!BoundNodeFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  BoundElemFile; BoundElemFile.open(BoundElemFilePath, std::ios::in | std::ios::binary); if (!BoundElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  ElemGeoIdsFile; ElemGeoIdsFile.open(ElemGeoIdsFilePath, std::ios::in | std::ios::binary); if (!ElemGeoIdsFile.is_open()) throw_line("ERROR, can't open input data file");
    
    //-------------------------
    // READ MESH FROM NODE FILE
    //-------------------------
    std::string line;
    std::istringstream iss;
    STREAM::getValue(NodeFile, line, iss, (*physics).nNodes);
    int nNodes = (*physics).nNodes;
    int dim = (*physics).dim;
    (*physics).nDof = nNodes;
    // choose number of thread dependding on the size of the problem
    unsigned int maxNumThreads = std::thread::hardware_concurrency();
    if ((1e5 <= (*physics).nDof) && ((*physics).nDof < 5e5))
    {
        PARALLEL::nThread = 8;
    }
    else if ((5e5 <= (*physics).nDof) && ((*physics).nDof < 1e6))
    {
        PARALLEL::nThread = 14;
    }
    else if (1e6 <= (*physics).nDof)
    {
        PARALLEL::nThread = maxNumThreads;
    }
    
    (*physics).coord.initialize(nNodes, dim);
    //---
    for(int i = 0; i < nNodes; i++)
    {
        prec* p = (*physics).coord[i];
        STREAM::getRowVector(NodeFile, line, iss, p, dim);
    }
    NodeFile.close();

    //-----------------------------
    // READ ELEMENTS FROM ELEM FILE
    //-----------------------------
    //---
    // pressure elements
    //---
    STREAM::getValue(ElemFile, line, iss, (*physics).nElem);
    int nElem = (*physics).nElem;
    int NodesxElem = dim+1;
    (*physics).elem.initialize(nElem, NodesxElem);

    //---
    for(int i = 0; i < nElem; i++)
    {
        int* p = (*physics).elem[i];
        STREAM::getRowVector(ElemFile, line, iss, p, NodesxElem);
    }
    ElemFile.close();
    
    //--------------------------------
    // READ BOUND INFO FROM BOUND FILE
    //--------------------------------
    std::getline(BoundElemFile, line);
    iss.str(line);
    iss >> nBoundElems; iss >> nBound;
    (*physics).nBounds = nBound;
    iss.clear();
    int NodesxBoundElem = dim;

    //----- alloc neuBoundIdElems ---
    neuBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
    neuBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
    neuBoundIdElems      = std::shared_ptr<int**[]> (new int**[nNeuBound]);
    nNeuBoundIdElems     = std::shared_ptr<int[]> (new int[nNeuBound]);
    //---
    int count = 0; int countEl = 0;
    int countNeu = 0;
    
    //----- alloc neuTimeBoundIdElems ---
    neuTimeBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
    neuTimeBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
    neuTimeBoundIdElems      = std::shared_ptr<int**[]> (new int**[nNeuBound]);
    nNeuTimeBoundIdElems     = std::shared_ptr<int[]> (new int[nNeuBound]);
    //---
    int countTime = 0; int countElTime = 0;
    int countNeuTime = 0;

    //----- alloc symmBoundIdElems ---
    symmBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
    symmBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
    symmBoundIdElems      = std::shared_ptr<int**[]> (new int**[nSymmBound]);
    nSymmBoundIdElems     = std::shared_ptr<int[]> (new int[nSymmBound]);
    
    //---
    int countS = 0; int countElSymm = 0;
    int countSymm = 0;
    
    for (int i = 0; i < nBound; i++)
    {
        int currID;
        int currNEl;
        std::getline(BoundElemFile, line);
        iss.str(line);
        iss >> currID; iss >> currNEl;
        iss.clear();
        if (neuTimeBound.hasIn(currID))
        {
            neuTimeBoundIdElems[countNeuTime] = &neuTimeBoundIdElems_buff[countElTime];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &neuTimeBoundIdElems_buff_buff[countTime];
                neuTimeBoundIdElems_buff[countElTime] = p;
                STREAM::getRowVector(BoundElemFile, line, iss, p, NodesxBoundElem);
                countTime += NodesxBoundElem;
                countElTime++;
            }
            nNeuTimeBoundIdElems[countNeuTime] = currNEl;
            countNeuTime++;
            std::getline(BoundElemFile, line);
        } 
        else if (neuBound.hasIn(currID))
        {
            neuBoundIdElems[countNeu] = &neuBoundIdElems_buff[countEl];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &neuBoundIdElems_buff_buff[count];
                neuBoundIdElems_buff[countEl] = p;

                STREAM::getRowVector(BoundElemFile, line, iss, p, NodesxBoundElem);
                count += NodesxBoundElem;
                countEl++;
            }
            nNeuBoundIdElems[countNeu] = currNEl;
            countNeu++;
            std::getline(BoundElemFile, line);
        } 
        else if (symmBound.hasIn(currID))
        {
            symmBoundIdElems[countSymm] = &symmBoundIdElems_buff[countElSymm];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &symmBoundIdElems_buff_buff[countS];
                symmBoundIdElems_buff[countElSymm] = p;
                STREAM::getRowVector(BoundElemFile, line, iss, p, NodesxBoundElem);
                countS += NodesxBoundElem;
                countElSymm++;
            }
            nSymmBoundIdElems[countSymm] = currNEl;
            countSymm++;
            std::getline(BoundElemFile, line);
        } 
        else
        {
            STREAM::getLines(BoundElemFile, line, currNEl+1);
            continue;
        }
    }
    // neuBoundIdElems_buff = (int**) realloc(neuBoundIdElems_buff, countEl*sizeof(int*));
    // neuBoundIdElems_buff_buff = (int*) realloc(neuBoundIdElems_buff_buff, count*sizeof(int));
    BoundElemFile.close();
    //----- alloc boundIdNodes ---
    std::getline(BoundNodeFile, line);
    iss.str(line);
    iss >> nBoundNodes;
    iss.clear();

    boundIdNodes_buff = std::shared_ptr<int[]> (new int[nBoundNodes]);
    boundIdNodes      = std::shared_ptr<int*[]> (new int*[nBound]); 
    nBoundIdNodes     = std::shared_ptr<int[]> (new int[nBound]);
    int countNodes = 0;
    for (int i = 0; i < nBound; i++)
    {
        int currID;
        int currNNodes;
        std::getline(BoundNodeFile, line);
        iss.str(line);
        iss >> currID; iss >> currNNodes;
        iss.clear();
        boundIdNodes[i] = &boundIdNodes_buff[countNodes];
        for (int j = 0; j < currNNodes; j++)
        {
            STREAM::getValue(BoundNodeFile, line, iss, boundIdNodes_buff[countNodes]);
            countNodes++;
        }
        std::getline(BoundNodeFile, line);
        nBoundIdNodes[i] = currNNodes;
    }
    BoundNodeFile.close();

    //--------------------------------------------
    // RECOVER PRESSURE BOUNDARY ELEMS FOR THE SYMMETRY CONDITIONS
    //--------------------------------------------
    BoundElemFile.open(BoundElemFilePath, std::ios::in | std::ios::binary); if (!BoundElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::getline(BoundElemFile, line);
    NodesxBoundElem = dim;
    //----- alloc ---
    symmPressureBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
    symmPressureBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
    symmPressureBoundIdElems      = std::shared_ptr<int**[]> (new int**[nSymmBound]);

    dirPressureBoundIdElems_buff_buff = std::shared_ptr<int[]> (new int[nBoundElems * NodesxBoundElem]); 
    dirPressureBoundIdElems_buff = std::shared_ptr<int*[]> (new int*[nBoundElems]);
    dirPressureBoundIdElems      = std::shared_ptr<int**[]> (new int**[nDirBound]);
    //---
    countS = 0; countElSymm = 0;
    countSymm = 0;

    int countD = 0; int countElDir = 0;
    int countDir = 0; dirBoundPressureNEl.initialize(nDirBound);
    
    for (int i = 0; i < nBound; i++)
    {
        int currID;
        int currNEl;
        std::getline(BoundElemFile, line);
        iss.str(line);
        iss >> currID; iss >> currNEl;
        iss.clear();
        if (symmBound.hasIn(currID))
        {
            symmPressureBoundIdElems[countSymm] = &symmPressureBoundIdElems_buff[countElSymm];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &symmPressureBoundIdElems_buff_buff[countS];
                symmPressureBoundIdElems_buff[countElSymm] = p;
                STREAM::getRowVector(BoundElemFile, line, iss, p, NodesxBoundElem);
                countS += NodesxBoundElem;
                countElSymm++;
            }
            countSymm++;
            std::getline(BoundElemFile, line);
        } 
        else if (dirBound.hasIn(currID))
        {
            dirPressureBoundIdElems[countDir] = &dirPressureBoundIdElems_buff[countElDir];
            dirBoundPressureNEl[countDir] = currNEl;
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &dirPressureBoundIdElems_buff_buff[countD];
                dirPressureBoundIdElems_buff[countElDir] = p;
                STREAM::getRowVector(BoundElemFile, line, iss, p, NodesxBoundElem);
                countD += NodesxBoundElem;
                countElDir++;
            }
            countDir++;
            std::getline(BoundElemFile, line);
        }
        else
        {
            STREAM::getLines(BoundElemFile, line, currNEl+1);
            continue;
        }
    }
    BoundElemFile.close();

    int n_geo_entities_ids;
    STREAM::getValue(ElemGeoIdsFile, line, iss, n_geo_entities_ids);
    if (n_geo_entities_ids != nElem) throw_line("ERROR: nElem and n_geo_entities_ids are different.\n");
    (*physics).elem_geo_entities_ids.initialize(nElem);
    STREAM::getColVector(ElemGeoIdsFile, line, iss, (*physics).elem_geo_entities_ids, nElem); 
    ElemGeoIdsFile.close();
    (*physics).max_geo_entity_id = (*physics).elem_geo_entities_ids.max() + 1;

    //--------------------------------------------
    // ALLOC EVERYTHING USEFUL FOR THE PROBLEM
    //--------------------------------------------
    nDof = nNodes;
    rhs_statForcing.setZeros(nDof); 
    rhs_timeForcing.initialize(nDof); 
    rhs.initialize(nDof);
}

// SET NODAL PERMEABILITIES
void PROBLEM_DARCY::set_permeabilities()
{
    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    int nNodes = (*physics).nNodes;

    // set discrete permeabilities
    discrete_permeabilities.setZeros(nNodes);
    VECTOR_INT nodal_priorities;
    nodal_priorities.setZeros(nNodes);
    for (int iel = 0; iel < nElem; iel++)
    {
        int* tempElem = (*physics).elem[iel];
        int idom = (*physics).elem_geo_entities_ids[iel];
        prec el_permeability = domains_permeability[idom];
        int permability_priority = domains_permeability_priority.getFirstOccurrence(idom, true);
        for (int inod = 0; inod < (dim+1); inod++)
        {
            int iglob = tempElem[inod];
            if (permability_priority >= nodal_priorities[iglob])
            {
                nodal_priorities[iglob] = permability_priority;
                discrete_permeabilities[iglob] = el_permeability;
            }            
        }
    }

    // set smooth permeabilities
    smooth_permeabilities.setZeros(nNodes);
    VECTOR nodal_weights;
    nodal_weights.setZeros(nNodes);
    for (int iel = 0; iel < nElem; iel++)
    {
        int* tempElem = (*physics).elem[iel];
        int idom = (*physics).elem_geo_entities_ids[iel];
        prec el_permeability = domains_permeability[idom];
        prec temp_vol = (*physics).Volume[iel];
        for (int inod = 0; inod < (dim+1); inod++)
        {
            int iglob = tempElem[inod];
            nodal_weights[iglob] += temp_vol;
            smooth_permeabilities[iglob] += el_permeability * temp_vol;
        }
    }
    smooth_permeabilities /= nodal_weights;
}

void PROBLEM_DARCY::setBC()
{
    if (completeLog < 2) std::cout << "\n----------\n--| SET BC |--\n----------\n";
    // check BC conditions
    for (int i = 0; i < nNoFluxBound; i++)
    {
        int bound = noFluxBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: noFlux geometry index out of bounds.\n");
    }
    for (int i = 0; i < nDirTimeBound; i++)
    {
        int bound = dirTimeBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Temporal Dirichlet geometry index out of bounds.\n");
    }
    for (int i = 0; i < nNeuTimeBound; i++)
    {
        int bound = neuTimeBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Temporal Neumann geometry index out of bounds.\n");
    }
    for (int i = 0; i < nDirBound; i++)
    {
        int bound = dirBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Static Dirichlet geometry index out of bounds.\n");
    }
    for (int i = 0; i < nNeuBound; i++)
    {
        int bound = neuBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Static Neumann geometry index out of bounds.\n");
    }
    for (int i = 0; i < nInnerBound; i++)
    {
        int bound = innerBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Inner Bound geometry index out of bounds.\n");
    }
    for (int i = 0; i < nSymmBound; i++)
    {
        int bound = symmBound[i];
        if (bound > nBound-1) throw_line("USAGE ERROR: Symmetry geometry index out of bounds.\n");
    }
    
    //-----------------------------------------------
    // CHECK THAT ALL BOUNDARY IDs HAVE BEEN DEFINED
    //-----------------------------------------------
    if (flagNoFlux == 0)
    {
        for (int ibound = 0; ibound < nBound; ibound++)
        {
            if (noFluxBound.hasIn(ibound)) continue;
            if (dirBound.hasIn(ibound)) continue;
            if (neuBound.hasIn(ibound)) continue;
            if (dirTimeBound.hasIn(ibound)) continue;
            if (neuTimeBound.hasIn(ibound)) continue;
            if (innerBound.hasIn(ibound)) continue;
            if (symmBound.hasIn(ibound)) continue;
            std::string errormsg = "ERROR: geometry index " + std::to_string(ibound) + " not defined in any type of BC";
            throw_line(errormsg);
        }
    }
    else if (flagNoFlux == 1)
    {
        nNoFluxBound = nBound - nDirBound - nNeuBound - nDirTimeBound - nNeuTimeBound - nInnerBound - nSymmBound;
        noFluxBound.setZeros(nNoFluxBound);
        int count = 0;
        for (int ibound = 0; ibound < nBound; ibound++)
        {
            if (dirBound.hasIn(ibound))
            {
                continue;
            } 
            if (neuBound.hasIn(ibound))
            {
                continue;
            } 
            if (dirTimeBound.hasIn(ibound))
            {
                continue;
            }
            if (neuTimeBound.hasIn(ibound))
            {
                continue;
            }
            if (innerBound.hasIn(ibound))
            {
                continue;
            }
            if (symmBound.hasIn(ibound))
            {
                continue;
            }
            noFluxBound[count] = ibound;
            count++;
        }
    }
    
    // initialize boundInfoMat 
    int nNodes = (*physics).nNodes;
    MATRIX_INT boundInfoMat(2, nNodes);

    // set bound nodes entries to -1
    for (int ibound = 0; ibound < nBound; ibound++)
    {
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            boundInfoMat[0][tempNode] = 0;
        }
    }

    //--- ITERATE ON THE TIME DIR BOUNDARIES ----
    dirTimeIdCount.initialize(nDirTimeBound);
    int sum = 0; 
    for (int idir = 0; idir < nDirTimeBound; idir++)
    {
        int ibound = dirTimeBound[idir];
        sum += nBoundIdNodes[ibound];
    }
    dirTimeNod.initialize(sum); //initialize dirNod with max possible dim
    nTimeDir = 0;
    for (int idir = 0; idir < nDirTimeBound; idir++)
    {
        int ibound = dirTimeBound[idir];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = 4;
                boundInfoMat[1][tempNode] = nTimeDir;
                dirTimeNod[nTimeDir] = tempNode;
                nTimeDir++;
            }
        }
        dirTimeIdCount[idir] = nTimeDir;
    }
    //dirTimeNod.shrink(nTimeDir);
    dirTimeNod.length = nTimeDir;
    int dim = (*physics).dim;
    dirTimeVal.initialize(nTimeDir, 1);

        //--- ITERATE ON THE STATIC DIR BOUNDARIES ----
    dirIdCount.initialize(nDirBound);
    sum = 0; 
    for (int idir = 0; idir < nDirBound; idir++)
    {
        int ibound = dirBound[idir];
        sum += nBoundIdNodes[ibound];
    }
    dirNod.initialize(sum); //initialize dirNod with max possible dim
    nDir = 0;
    for (int idir = 0; idir < nDirBound; idir++)
    {
        int ibound = dirBound[idir];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = 3;
                boundInfoMat[1][tempNode] = nDir;
                dirNod[nDir] = tempNode;
                nDir++;
            }
        }
        // printf("nBoundNodesDir of bound %d: %d\n", ibound, nBoundIdNodes[ibound]);
        dirIdCount[idir] = nDir;
    }
    //dirNod.shrink(nDir);
    dirNod.length = nDir;
    dirVal.initialize(nDir, 1);


    //--- ITERATE ON THE TIME NEU BOUNDARIES ----
    sum = 0; 
    for (int ineu = 0; ineu < nNeuTimeBound; ineu++)
    {
        sum += nNeuTimeBoundIdElems[ineu];
    }
    sum*=dim; // at most dim*#of bound elems of type Time Neu 
    neuTimeIdCount.initialize(nNeuTimeBound);
    neuTimeNod.initialize(sum); //initialize neuTimeNod with max possible dim
    neuTimeBaseFlux.initialize(sum); 
    neuTimeNormal.initialize(sum, dim);
    nTimeNeu = 0;
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    for (int ineu = 0; ineu < nNeuTimeBound; ineu++)
    {
        
        int** tempNeuBoundElems = neuTimeBoundIdElems[ineu];
        for (int iel = 0; iel < nNeuTimeBoundIdElems[ineu]; iel++)
        {
            MATRIX matCoord(dim,dim); 
            VECTOR normal(dim);
            int* tempElem = tempNeuBoundElems[iel];

            for (int inode = 0; inode < dim; inode++)
            {
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord[tempElem[inode]][icomp]; 
            }
            prec tempBaseFlux = getArea(matCoord, dim)/dim;

            for (int inode = 0; inode < dim; inode++)
            {
                int tempNode = tempElem[inode];

                if (boundInfoMat[0][tempNode] == 0)
                {
                    boundInfoMat[0][tempNode] = 2;
                } 

                if (boundInfoMat[0][tempNode] == 2)
                {
                    neuTimeNod[nTimeNeu] = tempNode;
                    neuTimeBaseFlux[nTimeNeu] = tempBaseFlux;
                    nTimeNeu++;
                }
            };
        }
        neuTimeIdCount[ineu] = nTimeNeu;
    }
    neuTimeNod.shrink(nTimeNeu); neuTimeBaseFlux.shrink(nTimeNeu);
    neuTimeNod.length = nTimeNeu; neuTimeBaseFlux.length = nTimeNeu;
    neuTimeVal.zeros(nTimeNeu, 1);

    //--- ITERATE ON THE STATIC NEU BOUNDARIES ----
    sum = 0; 
    for (int ineu = 0; ineu < nNeuBound; ineu++)
    {
        int ibound = neuBound[ineu];
        sum += nBoundIdNodes[ibound];
    }
    neuNod.initialize(sum); //initialize neuNod with max possible dim
    nNeu = 0;
    for (int ineu = 0; ineu < nNeuBound; ineu++)
    {
        int ibound = neuBound[ineu];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = 1;
                boundInfoMat[1][tempNode] = nNeu;
                neuNod[nNeu] = tempNode;
                nNeu++;
            }
        }
    }
    neuNod.shrink(nNeu);
    neuNod.length = nNeu;
    neuVal.zeros(nNeu, 1);

    //--- ITERATE ON THE SYMMETRY BOUNDARIES ----
    sum = 0; 
    for (int isymm = 0; isymm < nSymmBound; isymm++)
    {
        int ibound = symmBound[isymm];
        sum += nBoundIdNodes[ibound];
    }
    symmNod.initialize(sum); //initialize neuSymm with max possible dim
    nSymm = 0;
    for (int isymm = 0; isymm < nSymmBound; isymm++)
    {
        int ibound = symmBound[isymm];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = -1;
                boundInfoMat[1][tempNode] = nSymm;
                symmNod[nSymm] = tempNode;
                nSymm++;
            }
        }
    }
    symmNod.shrink(nSymm);
    symmRealIglob.initialize(nSymm);
    symmNod.length = nSymm;
    symmVal.initialize(nSymm, dim);

    sum = 0; 
    for (int inoFlux = 0; inoFlux < nNoFluxBound; inoFlux++)
    {
        int ibound = noFluxBound[inoFlux];
        sum += nBoundIdNodes[ibound];
    }
    noFluxNod.initialize(sum); //initialize noFluxNod with max possible dim
    nNoFlux = 0;
    noFluxIdCount.initialize(nNoFluxBound);
    for (int inoFlux = 0; inoFlux < nNoFluxBound; inoFlux++)
    {
        int ibound = noFluxBound[inoFlux];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = -1;
                noFluxNod[nNoFlux] = tempNode;
                nNoFlux++;
            }
        }
        noFluxIdCount[inoFlux] = nNoFlux;
    }
    noFluxNod.shrink(nNoFlux);
    noFluxNod.length = nNoFlux;

    setStatDirBC();
    setStatNeuBC(boundInfoMat);
    setSymmBC(boundInfoMat);
}
//----------------------------------------------------------------------------
// SET BOUNDARY CONDITIONS
//-----------------------------------------------------------------------------
void PROBLEM_DARCY::setTimeDirBC(prec trialTime) //priority 4
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET TIME DIR BC |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    int old = 0;
    for (int iDirBound = 0; iDirBound < nDirTimeBound; iDirBound++)
    {
        int endPos = dirTimeIdCount[iDirBound];
        int tempNNodes = endPos - old;
        MATRIX tempCoord(tempNNodes, dim);
        int tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            int iglob = dirTimeNod[idir];
            prec* tempPoint = tempCoord[tempCount];
            prec* point  = coord[iglob];
            for (int icomp = 0; icomp < dim; icomp++) tempPoint[icomp] = point[icomp];
            tempCount++;
        }

        std::vector<VECTOR> tempTimeDirVal = dirTimeFunc[iDirBound](tempCoord, trialTime);

        tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            dirTimeVal[idir][0] = tempTimeDirVal[0][tempCount];
            tempCount++;
        }
        old = dirTimeIdCount[iDirBound];
    }
}
//---
void PROBLEM_DARCY::setStatDirBC() //priority 3
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| DARCY SET STAT DIR BC  |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    int old = 0;
    for (int iDirBound = 0; iDirBound < nDirBound; iDirBound++)
    {
        int endPos = dirIdCount[iDirBound];
        int tempNNodes = endPos - old;
        MATRIX tempCoord(tempNNodes, dim);
        int tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            int iglob = dirNod[idir];
            prec* tempPoint = tempCoord[tempCount];
            prec* point  = coord[iglob];
            for (int icomp = 0; icomp < dim; icomp++) tempPoint[icomp] = point[icomp];
            tempCount++;
        }

        std::vector<VECTOR> tempDirVal(1);
        FUNCTION_PARSER::evalCoord(dirFunc[iDirBound][0], tempCoord, 0.0, tempDirVal[0]);

        tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            prec* tempDir = dirVal[idir];
            tempDir[0] = tempDirVal[0][tempCount];
            tempCount++;
        }
        old = endPos;
    }
}
//---
void PROBLEM_DARCY::setTimeNeuBC(prec trialTime) // priority 2
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| DARCY SET TIME NEU BC |--\n----------\n";
    int old = 0;
    int neuCount = 0;
    for (int iNeuBound = 0; iNeuBound < nNeuTimeBound; iNeuBound++)
    {
        prec meanFlux = *neuTimeMeanFlux[iNeuBound](0, 0, 0, trialTime);
        for (int ineu = old; ineu < neuTimeIdCount[iNeuBound]; ineu++)
        {
            neuTimeVal[ineu][0] = neuTimeBaseFlux[neuCount] * meanFlux;
            neuCount++;
        }
        old = neuTimeIdCount[iNeuBound];
    }
}
//---
void PROBLEM_DARCY::setStatNeuBC(MATRIX_INT& boundInfoMat) // priority 1
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| DARCY SET STAT NEU BC |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    int* priority = boundInfoMat[0];
    int* neuNodPos = boundInfoMat[1];
    MATRIX matCoord(dim,dim); 
    VECTOR normal(dim);

    for (int ibound = 0; ibound < nNeuBound; ibound++)
    {
        prec meanFluxByDim = neuMeanFlux[ibound]/dim;
        //-------------------------------------
        for (int iel = 0; iel < nNeuBoundIdElems[ibound]; iel++)
        {
            int* el = neuBoundIdElems[ibound][iel];
            for (int inode = 0; inode < dim; inode++)
            {
                int iglob = el[inode];
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord[iglob][icomp]; 
            }
            prec Area = getArea(matCoord, dim);
            //---
            prec tempFactor = Area * meanFluxByDim;
            for (int iloc = 0; iloc < dim; iloc++)
            {
                int iglob = el[iloc];
                // std::cout << "iglob: " << iglob << "\n";
                if (priority[iglob] != 1) 
                {
                    continue;
                }
                int tempNeuNod = neuNodPos[iglob];
                neuVal[tempNeuNod][0] += tempFactor;
            }    
        }
    }
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
// SET SYMMETRY BC
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
void PROBLEM_DARCY::setSymmBC(MATRIX_INT& boundInfoMat) //priority -1
{
    // Symm BC are equivalent to NoFlux BC for this problem.
    //  Yhe NoFlux BC are equivalent to DoNothing. So nothing has to be done in this case.
}

//-----------------------------------------------------------------------
// FORCING 
//-----------------------------------------------------------------------
void PROBLEM_DARCY::applyStatForcing()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| DARCY APPLY STATIC FORCING |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes = (*physics).nNodes;
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
    int nElem = (*physics).nElem;

    std::vector<VECTOR> forcingAtNodes(1);
    forcingAtNodes[0].initialize(nNodes);
    FUNCTION_PARSER::evalCoord(statForcing[0], (*physics).coord, 0, forcingAtNodes[0]);
    //----------
    for (int iel = 0; iel < nElem; iel++)
    {
        int* tempElem = elem[iel];
        prec element_factor = Volume[iel] / ((dim+1)*(dim+2));
        prec tempForcingSum = 0;
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = tempElem[iloc];
            tempForcingSum += forcingAtNodes[0][iglob];
        }
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = tempElem[iloc];
            prec forcingFactor = (forcingAtNodes[0][iglob] + tempForcingSum);
            rhs_statForcing[iglob] += forcingFactor * element_factor;
        }
    }        
    //---------------------------------------
}
//---

void PROBLEM_DARCY::updateTimeForcing()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS UPDATE TIME FORCING |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes = (*physics).nNodes;
    int nElem = (*physics).nElem;
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
    prec evalTime = time + deltaT;
    rhs_timeForcing.resetZeros();
    std::vector<VECTOR> forcingAtNodes(1);

    forcingAtNodes[0].initialize(nNodes);
    FUNCTION_PARSER::evalCoord(timeForcing[0], (*physics).coord, evalTime, forcingAtNodes[0]);
    //----------
    
    for (int iel = 0; iel < nElem; iel++)
    {
        int* tempElem = elem[iel];
        prec element_factor = Volume[iel] / ((dim+1)*(dim+2));
        prec tempForcingSum = 0;
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            
            int iglob = tempElem[iloc];                
            tempForcingSum += forcingAtNodes[0][iglob];
            
        }
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = tempElem[iloc];
            prec forcingFactor = (forcingAtNodes[0][iglob] + tempForcingSum);
            rhs_timeForcing[iglob] += forcingFactor * element_factor;
        }
    }        
}
//------------------------------------------------------------------------------
// IMPOSE BOUNDARY CONDITIONS
//------------------------------------------------------------------------------
void PROBLEM_DARCY::imposeBC(CSRMAT &SYSMAT_final, VECTOR &rhs)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| IMPOSE BC |--\n----------\n";
    
    imposeStaticBC(SYSMAT_final, rhs);
    if ((*physics).isStationary == 0) imposeTimeBC(SYSMAT_final, rhs);
}

//------------------------------
// IMPOSE STATIC BC
//------------------------------
void PROBLEM_DARCY::imposeStaticBC(CSRMAT &SYSMAT_Darcy, VECTOR &rhs)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| IMPOSE STATIC BC  |--\n----------\n";
    int flagBC = (*physics).flagBC;

    // DIRICHLET BC
    switch (flagBC)
    {
        case 0:
        {
            prec penalty = 1e13;
            // noFlux: DO NOTHING
            // symm: DO NOTHING
            // static dir
            for (int idir = 0; idir < nDir; idir++)
            {
                int iglob = dirNod[idir];
                SYSMAT_Darcy(iglob,iglob) = penalty;
                rhs[iglob] = penalty * dirVal[idir][0];
            }
            break;
        }
        case 1:
        {
            
            if (!statBCAlreadyApplied)
            {
                // noFlux: DO NOTHING
                // symm: DO NOTHING
                
                // static dir
                // set to zero row & col entries
                int countPosToZero = 0; int countIglobToZero = 0; int countPosToOne = 0;
                std::shared_ptr<prec[]> coefMat = SYSMAT_Darcy.coef;
                dirCountToZero[0] = 0;
                for (int idir = 0; idir < nDir; idir++)
                {
                    prec tempDirVal = dirVal[idir][0];
                    int tempDirNod  = dirNod[idir];
                    rhs[tempDirNod] = tempDirVal;
                    VECTOR_INT rowInColToZero = SYSMAT_Darcy.setColInRowToZero(tempDirNod, dirPosToZero, countPosToZero, dirPosToOne, countPosToOne);
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int iglob = rowInColToZero[i];
                        int pos = SYSMAT_Darcy.getPos(iglob, tempDirNod);
                        prec tempVal = coefMat[pos];
                        rhs[iglob] -= tempVal * tempDirVal;
                        coefMat[pos] = 0.0;
                        dirIglobToZero[countIglobToZero] = iglob; countIglobToZero++;
                        dirPosToZero[countPosToZero] = pos; countPosToZero++;

                        int tempDirCountId = idir + 1;
                        dirCountToZero[tempDirCountId] = countPosToZero;
                    }
                }
                statBCAlreadyApplied = true;
            } 
            else
            {
                //--- Dirichlet
                int nToOne = dirPosToOne.length; int countIglobToZero = 0;
                std::shared_ptr<prec[]> coefMat = SYSMAT_Darcy.coef;
                for (int idir = 0; idir < nDir; idir++)
                {
                    int tempDirNod = dirNod[idir];
                    prec tempVal = dirVal[idir][0]; 
                    rhs[tempDirNod] = tempVal;
                    int tempDirCountId = idir;
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
                //---
                for (int i = 0; i < nToOne; i++)
                {
                    int pos = dirPosToOne[i];
                    coefMat[pos] = 1;
                }
            }
            //---
            break;
        }
    }

    // NEUMANN BC
    //static
        for (int ineu = 0; ineu < nNeu; ineu++)
        {
            int iglob = neuNod[ineu];
            rhs[iglob] -= neuVal[ineu][0];
        }
}

//-----------------------------
// IMPOSE TIME DEPENDENT BC
//-----------------------------
void PROBLEM_DARCY::imposeTimeBC(CSRMAT &SYSMAT_Darcy, VECTOR &rhs)
{
    static bool timeBCAlreadyApplied = false;
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS IMPOSE TIME BC  |--\n----------\n";
    int flagBC = (*physics).flagBC;
    // DIRICHLET BC
    switch (flagBC)
    {
        case 0:
        {
            prec penalty = 1e13;
            // time dep dir
            for (int idir = 0; idir < nTimeDir; idir++)
            {
                int iglob = dirTimeNod[idir];  
                SYSMAT_Darcy(iglob,iglob,penalty);
                rhs[iglob] = penalty * dirTimeVal[idir][0];
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
                std::shared_ptr<prec[]> coefMat = SYSMAT_Darcy.coef;
                dirTimeCountToZero[0] = 0;
                for (int idir = 0; idir < nTimeDir; idir++)
                {
                    prec tempDirVal = dirTimeVal[idir][0];
                    int tempDirTimeNod  = dirTimeNod[idir];
                    rhs[tempDirTimeNod] = tempDirVal;
                    VECTOR_INT rowInColToZero = SYSMAT_Darcy.setColInRowToZero(tempDirTimeNod, dirTimePosToZero, countPosToZero, dirTimePosToOne, countPosToOne); //since the matrix has symmetric non zero entries we use the                                                                    //non 
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int iglob = rowInColToZero[i];
                        int pos = SYSMAT_Darcy.getPos(iglob, tempDirTimeNod);
                        prec tempVal = coefMat[pos];
                        rhs[iglob] -= tempVal * tempDirVal;
                        coefMat[pos] = 0.0;
                        dirTimeIglobToZero[countIglobToZero] = iglob; countIglobToZero++;
                        dirTimePosToZero[countPosToZero] = pos; countPosToZero++;

                        int tempDirTimeCountId = idir + 1;
                        dirTimeCountToZero[tempDirTimeCountId] = countPosToZero;
                    }
                }
                timeBCAlreadyApplied = true;
            } 
            else
            {
                int nToOne = dirTimePosToOne.length;
                std::shared_ptr<prec[]> coefMat = SYSMAT_Darcy.coef; int countIglobToZero = 0;
                //---
                for (int idir = 0; idir < nTimeDir; idir++)
                {
                    int tempDirTimeNod = dirTimeNod[idir];
                    prec tempVal = dirTimeVal[idir][0];
                    rhs[tempDirTimeNod] = tempVal;

                    int tempDirTimeCountId = idir;
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
    for (int ineu = 0; ineu < nTimeNeu; ineu++)
    {
        int iglob = neuTimeNod[ineu];
        rhs[iglob] -= neuTimeVal[ineu][0];
    }
}

//-------------------------------------------------
// UPDATE BC
//-------------------------------------------------
void PROBLEM_DARCY::updateBC(prec trialTime)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS UPDATE BC  |--\n----------\n";
    setTimeDirBC(trialTime);
    setTimeNeuBC(trialTime);
}

//-------------------------------------------------
// UPDATE SYSMAT_Darcy 
//-------------------------------------------------
void PROBLEM_DARCY::updateSYSMAT()
{  
    // add Time Derivative Matrix to SYSMAT
    prec factorP = 1 / deltaT;
    addToSysmat(M, factorP); 
}
//---

//-------------------------------------------------
// UPDATE RHS
//-------------------------------------------------
void PROBLEM_DARCY::updateRHS()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| UPDATE RHS |--\n----------\n";
    int nNodes = (*physics).nNodes;
    rhs.completeReset();
    rhs = rhs_statForcing;

    if (flagForcing == 1)
    {
        updateTimeForcing();
        rhs += rhs_timeForcing;
    }

    VECTOR last_sol_contribution(nNodes);
    last_sol_contribution = (M*lastSol) / deltaT;
    rhs += last_sol_contribution;
}
//------------------------------------------------------------------------------------
void PROBLEM_DARCY::setInitCond()
{
    // define init Sol 

    lastSol.setZeros(nDof);


    // // impose stat dir

    // for (int idir = 0; idir < nDir; idir++)
    // {
    //     int iglob = dirNod[idir];

    //     prec* tempVal = dirVal[idir];
    //     for (int icomp = 0; icomp < dim; icomp++)
    //     {
    //         int globNod = iglob+icomp*nNodes_v;
    //         lastSol[globNod] = tempVal[icomp];
    //     }
    // }
}

//---------------------------------------------------
// PREPARE LINEAR SYSTEM (ASSEMBLY OF THE MATRICES)
//---------------------------------------------------
void PROBLEM_DARCY::prepareSolver()
{
    assemble();

    applyStatForcing();

    // prepare BC update in case of lifting functions imposition
    int dim = (*physics).dim;
    int nElem = (*physics).nElem;
    std::shared_ptr<int*[]> elem = (*physics).elem.PP;
    if ((*physics).flagBC == 1)
    {
        std::shared_ptr<int[]> iat = SYSMAT_base.iat;
        int countPosToZero = 0;
        //--- noFlux
        for (int inoFlux = 0; inoFlux < nNoFlux; inoFlux++)
        {
            int iglob = noFluxNod[inoFlux];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nNoFlux;
        countPosToZero *= 2;
        noFluxPosToZero.initialize(countPosToZero);
        noFluxPosToOne.initialize(nNoFlux);
        //--- static
        countPosToZero = 0;
        for (int idir = 0; idir < nDir; idir++)
        {
            int iglob = dirNod[idir];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nDir;
        dirIglobToZero.initialize(countPosToZero);
        countPosToZero *= 2;
        dirPosToZero.initialize(countPosToZero);
        dirPosToOne.initialize(nDir);
        dirCountToZero.initialize(nDir + 1);
        //--- time-dependent
        countPosToZero = 0;
        for (int idir = 0; idir < nTimeDir; idir++)
        {
            int iglob = dirTimeNod[idir];
            int nnzInRow = iat[iglob+1] - iat[iglob];
            countPosToZero += nnzInRow;
        }
        countPosToZero -= nTimeDir;
        dirTimeIglobToZero.initialize(countPosToZero);
        countPosToZero *= 2;
        dirTimePosToZero.initialize(countPosToZero);
        dirTimePosToOne.initialize(nTimeDir);
        dirTimeCountToZero.initialize(nTimeDir + 1);

        //--------------------------------------------------
        // SYMMETRY and NO-FLUX BC
        //--------------------------------------------------
        // symm BC are equivalent to "do-nothing"
    }
    
    //-----------------------
    // FREE USELESS STUFF
    //-----------------------
    boundIdNodes_buff.reset(); boundIdNodes.reset(); 
    // // dirBound.dlt();
    // // noFluxBound.dlt(); 

    //----------------------------------------------------------
    // get max length of each element
    //----------------------------------------------------------
    (*physics).h.initialize(nElem);
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    std::vector<VECTOR> tempCoord(dim+1);
    for (int iel = 0; iel < nElem; iel++)
    {
        prec maxDist = 0;
        int* tempElem = elem[iel];
        //----------------------------------------
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = tempElem[iloc];
            prec* nodeCoord = coord[iglob];
            tempCoord[iloc].initialize(dim);
            for (int icomp = 0; icomp < dim; icomp++)
            {
                tempCoord[iloc][icomp] = nodeCoord[icomp];
            }
        }
        //----------------------------------------
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            for (int jloc = iloc+1; jloc < dim+1; jloc++)
            {
                VECTOR diff = tempCoord[iloc] - tempCoord[jloc];
                prec tempDist = diff.norm();
                if (tempDist > maxDist)
                {
                    maxDist = tempDist;
                }
            }
        }
        (*physics).h[iel] = maxDist;
    }

    //--------------------------------------
    // PRINT INIT COND
    //--------------------------------------
    int nNodes = (*physics).nNodes;
    std::string folderName = "Darcy_sol";
    VTKWriter.initializeForNS(name, dim, nNodes, nElem, folderName);

    setInitCond();
}

//------------------------------
// RESET PRINT
//------------------------------
void PROBLEM_DARCY::resetPrint(int iter)
{
    deltaT = deltaT0;
    if (save_solution == 1)
    {
        std::string newFolderName = "iter" + std::to_string(iter);
        VTKWriter.resetForNS(newFolderName);
        printRes = true;
    }
}

//-----------------------------------------------
// SOLVER
//-----------------------------------------------
//-------------------------
void PROBLEM_DARCY::StatSolver()
{
    time = 0;
    deltaT = (*physics).deltaT;
    globIter = 0;
    (*physics).Darcy_solution.complete_reset();

    oneStepSolver(1e-3, 10);
    (*physics).Darcy_solution.append_row(lastSol);
    print_sol_in_VTK((*physics).Darcy_solution);

    VTKWriter.closeTFile();
}
//----

void PROBLEM_DARCY::Solver()
{
    real_solution_times.resetZeros();
    real_solution_times.initialize(0);
    (*physics).Darcy_solution.complete_reset();

    // pos dal primo BC
    time = 0;
    real_solution_times.append(time);
    globIter = 0;
    prec t_end = (*physics).t_end;
    int nTimeSteps = t_end/deltaT;
    if (nTimeSteps > 1) lastSol.resetZeros();   // in the NS system there are the convergence iterations,
                                                // so in case of a single time step, starting from the lastSol obtained
                                                // can fasten the solution procedure starting from an initial condition similar to the solution.
                                                // Unfortunately, in case of more times steps it's not possible to exploit this "trick". 
                                                // Maybe it could be useful to recall the solution of the last optimization iteration 
                                                // and try to find every time step solution starting from a similar one obtained in the simulation before.
    printf("|-| IT: 0 \tTIME: %7.4" format " \tINITIAL CONDITION\n", 0.0);

    VECTOR init_cond((*physics).nDof);
    init_cond.setZeros(((*physics).nDof));
    (*physics).Darcy_solution.append_row(init_cond);

    while ((t_end-time) > (1e-15 * (nTimeSteps+1)))
    {
        oneStepSolver();
        (*physics).Darcy_solution.append_row(lastSol);
        if (time+deltaT > t_end) deltaT = t_end-time;
    }
    print_sol_in_VTK((*physics).NS_solution);
    VTKWriter.closeTFile();
}

//-------------------------
// ONE-STEP SOLVER
//-------------------------
void PROBLEM_DARCY::oneStepSolver(prec toll, int itMax)
{
    prec trialTime = time + deltaT;
    SYSMAT = SYSMAT_base; // initialize SYSMAT as SYSMAT_base in order to add updated matrices in correct blocks
    if ((*physics).isStationary == 0)
    {
        updateBC(trialTime);
        updateSYSMAT();
        updateRHS();
    }

    // if (completeLog) 
    if (completeLog < 2) std::cout << "\n----------------------------------------------------------------------------\n--| DARCY ONE-STEP SOLVER TIME " << time << "-> " << trialTime << "|--\n----------------------------------------------------------------------------\n";

    VECTOR rhs_final;
    CSRMAT SYSMAT_final;
    prec finalRes = 0.0;
    // prec scarto = 0.0;

    VECTOR oldSol;
    oldSol = lastSol;
    VECTOR sol(nDof);

    double startTime = 0;
    double endTime = 0;
    // double startSolverTime = 0;
    // double endSolverTime = 0;

    startTime = omp_get_wtime();
    rhs_final = rhs;
    SYSMAT_final = SYSMAT;

    imposeBC(SYSMAT_final, rhs_final);

    std::shared_ptr<prec []> solP = sol.P;

    // SOLVERLS::launchPardiso(SYSMAT_final, rhs_final.P, solP, 1); finalRes = toll-1;
    precondJacobiSolver(SYSMAT_final, rhs_final, sol);
    // precondSchurSolver(SYSMAT_final, rhs_final, sol);
    // SIMPLESolver(SYSMAT_final, rhs_final, sol, finalRes);
    // SCHURgmres(SYSMAT_final, rhs_final, sol, finalRes);
    // prec bNorm   = rhs_final.norm();
    // sol = SOLVERLS::gmres_p(SYSMAT_final, rhs_final.P, lastSol.P, bNorm, 1e-6, 500); 

    bool is_nan = false;
    for (int i = 0; i < nDof; i++) if(!(abs(sol[i]) < 1e16))
    {
        is_nan = true;
        break;
        // throw_line("THE ERROR IS HERE MAN\n");
    }
    if (is_nan) 
    {
        // sol = SOLVERLS::launchPardiso(SYSMAT_final, rhs_final.P, 0);

        // prec bNorm   = rhs_final.norm();
        // sol = SOLVERLS::gmres_p(SYSMAT_final, rhs_final.P, lastSol.P, bNorm, 1e-6, 500); 
    }
    // update last sol
    PARALLEL::copy(sol, lastSol);
    // scarto = evalErrRelV(oldSol.P, sol.P);
    //printf("TIME INNER |SIMPLE| SOLVER: %" format ", scarto: %e\n", endTime-startTime, scarto);
    // printf("scarto it %d è : %" format_e "\n", it, sol);

    PARALLEL::copy(lastSol, oldSol);

    endTime = omp_get_wtime();

    printf("|DARCY| IT: %d \tTIME: %7.4" format " \tSOLVER TIME: %" format ", residual: %e\n", (globIter+1), trialTime, endTime-startTime, finalRes);    
    time = trialTime;
    real_solution_times.append(time);
    prec t_end = (*physics).t_end;  
    if (completeLog < 2) printf("In progress... %%%4.2f\n", time/t_end*100);    
    globIter++;
}
//

//-------------------------------
// PRINT STEP SOLUTION IN VTK
//-------------------------------

void PROBLEM_DARCY::print_sol_in_VTK(MATRIX &requested_sol)
{
    for (int itime = 0; itime < (*physics).solution_times.length; itime++)
    {
        prec trial_time;
        if ((*physics).isStationary == 1) //stationary case: only one time = "infinity", so it has no meaning to report it.
        {
            trial_time = 0;
        }
        else // time dependent case: starting from t0.
        {
            trial_time = (*physics).solution_times[itime];
        }
         
        // GET PRESSURE
        VECTOR pressure = requested_sol.get_row(itime);

        std::cout << "\nhere\n";
        std::vector<VECTOR> pressure_gradient;
        (*physics).eval_gradient_p(pressure, pressure_gradient);
        // MATRIX velocity;
        // (*physics).convert_std_vector_of_vectors_into_matrix(pressure_gradient, velocity);
        // velocity /= (*physics).mu;
        // velocity *= smooth_permeabilities;
        // pause();
        VECTOR velocity_magnitude;
        (*physics).eval_gradient_norm(pressure_gradient, velocity_magnitude);
        velocity_magnitude /= (*physics).mu;
        velocity_magnitude *= smooth_permeabilities;
        
        VECTOR bound_ids(5);
        bound_ids[0] = 11; 
        bound_ids[1] = 12; 
        bound_ids[2] = 13; 
        bound_ids[3] = 14; 
        bound_ids[4] = 15; 
        eval_flux_sum_on_boundaries(bound_ids, velocity_magnitude);
        
        
        //-----------------
        if (printRes) VTKWriter.write((*physics).coord, (*physics).elem, trial_time, pressure, 1, "Pressure", discrete_permeabilities, 1, "zK_d", smooth_permeabilities, 1, "zK_s", velocity_magnitude, 1, "||Velocity||");
    }    
}

void PROBLEM_DARCY::eval_flux_sum_on_boundaries(VECTOR bound_ids, VECTOR &velocity)
{
    // N.B.: the velocity vector is supposed aligned with the normal to each element

    prec flux = 0;
    prec total_area = 0;

    int dim = (*physics).dim;
    for (int ibound = 0; ibound < bound_ids.length; ibound++)
    {
        int bound_id = bound_ids[ibound];
        MATRIX_INT bound_elems = (*physics).bound_elems[bound_id];
        for (int iel = 0; iel < bound_elems.nRow; iel++)
        {
            int* temp_elem = bound_elems[iel];
            MATRIX coord(dim, dim);
            for (int inode = 0; inode < dim; inode++)
            {
                int node = bound_elems[iel][inode];
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    coord[inode][jcomp] = (*physics).coord[node][jcomp];
                }
            }
            prec el_area = (*physics).get_surface(coord, dim);
            total_area += el_area;
            for (int inod = 0; inod < dim; inod++)
            {
                int iglob = temp_elem[inod];
                flux += velocity[iglob] * el_area / dim;
            }
        }
    }

    prec avg_velocity = flux / total_area;
    prec domain_length = 100;
    prec delta_P = 1000;
    prec eq_permeability = avg_velocity * domain_length * (*physics).mu / delta_P;

    std::cout << "\n   ----------------------------------------------------";
    std::cout << "\n   | DARCY PROBLEM EQUIVALENT PERMEABILITY |";
    std::cout << "\n   ---| REQUESTED FLUX: " << flux << "\n";
    std::cout << "   ---| TOTAL AREA: " << total_area << "\n";
    std::cout << "   ---| AVG. VELOCITY: " << avg_velocity << "\n";
    std::cout << "   ---| EQ. PERMEABILITY: " << eq_permeability << "\n";
    std::cout << "   ----------------------------------------------------\n";
    
}

// PRECONDITIONER !!!!!!!
void PROBLEM_DARCY::precondJacobiSolver(CSRMAT &SYS,VECTOR &rhs, VECTOR &sol)
{

    VECTOR diag(nDof); std::shared_ptr<prec[]> diagP = diag.P;

    //  PRECOND THE SYSTEM MATRIX
    SYS.precondNSDiag(diagP); diag.invert();
    CSRMAT jacSYS;
    PARALLEL::RprodByDiag(SYS, diag, jacSYS);

    VECTOR precondLastSol(nDof); PARALLEL::pointProd(precondLastSol, diag, precondLastSol); 
    prec bNorm = PARALLEL::norm(rhs);

    PARALLEL::gmres_p(jacSYS, rhs, sol, precondLastSol, bNorm, 1e-5, 500); 

    PARALLEL::pointProd(sol, diag, sol);
}

void PROBLEM_DARCY::checkImportParameters()
{
    std::cout << "\n+++++++++++++++| CHECK PHYSICS IMPORT-PARAMETERS |+++++++++++++++\n\n";

    // PROBLEM NAME & DIMENSION
    std::cout << " Problem Name: " << name << "\n" ;
    std::cout << " Dimension: " << (*physics).dim << "\n";

    // TIME PARAMETERS
    std::cout << "\n.-=| TIME PARAMETERS |=-. \n";
    std::cout << " Ending time:  " << (*physics).t_end << "\n" ;
    std::cout << " Default time step: " << deltaT << "\n" ;
    std::cout << " Minimal time step: " << (*physics).deltaT_min << "\n" ;

    // PHYSICS PARAMETERS
    std::cout << "\n.-=| PHYSICS PARAMETERS |=-. \n";
    std::cout << " rho: "<< (*physics).rho << "\n";
    std::cout << " mu: " << (*physics).mu << "\n" ;
    // PHYSICS PARAMETERS
    std::cout << "\n.-=| PERMEABILITY |=-. \n";
    std::cout << " n_domains: "<< n_domains << "\n";
    domains_permeability.printRowMatlab(" Domains Permeability");
    
    // FORCING
    std::cout << "\n.-=| FORCING |=-. \n";
    // forcing
    CUSTOM::printRowStd(statForcing, " STATIC Forcing: ");
    if (flagForcing == 1) CUSTOM::printRowStd(timeForcing, " TIME DEPENDENT Forcing: ");

    // BOUNDARY CONDITIONS
    std::cout << "\n.-=| BOUNDARY CONDITIONS |=-. \n";
    std::cout << " Flag BC: " << (*physics).flagBC << "\n";
    std::cout << " Flag No-Flux: " << flagNoFlux << "\n";
    // noFlux
    std::cout << "\n--| NO-FLUX:  ";
    std::cout << "\n # NO-FLUX Bound ID: " << nNoFluxBound << "\n";
    noFluxBound.printRowMatlab(" NO-FLUX Bound ID");
    // Inner 
    std::cout << "\n--| INNER:  ";
    std::cout << "\n # Inner Bound ID: " << nInnerBound << "\n";
    innerBound.printRowMatlab(" Inner Bound ID");
    // Symmetry
    std::cout << "\n--| SYMMETRY:  ";
    std::cout << "\n # Symmetry Bound ID: " << nSymmBound << "\n";
    symmBound.printRowMatlab(" Symmetry Bound ID");
    // Static Dirichlet
    std::cout << "\n--| STATIC DIRICHLET:  ";
    std::cout << "\n # STATIC Dirichlet Bound ID: " << nDirBound << "\n";
    dirBound.printRowMatlab(" STATIC Dirichlet Bound ID");
    for (int ibound = 0; ibound < nDirBound; ibound++)
    {
        std::string tempName = " " + std::to_string(ibound+1) + "th STATIC Dirichlet function: ";
        CUSTOM::printRowStd(dirFunc[ibound], tempName);
    }
    // Static Neumann
    std::cout << "\n--| STATIC NEUMANN:  ";
    std::cout << "\n # STATIC Neumann Bound ID: " << nNeuBound << "\n";
    neuBound.printRowMatlab(" STATIC Neumann Bound ID");
    neuMeanFlux.printRowMatlab(" STATIC Neumann flux");

    // Time Dependent Dirichlet
    std::cout << "\n--| TIME DEPENDENT DIRICHLET:  ";
    std::cout << "\n # TIME DEPENDENT Dirichlet Bound ID: " << nDirTimeBound << "\n";
    dirTimeBound.printRowMatlab(" TIME DEPENDENDT Dirichlet Bound ID");
    for (int i = 0; i < nDirTimeBound; i++)
    {
        std::string tempName = std::to_string(i+1) + "th TIME DEPENDENT Dirichlet function";
        dirTimeFunc[i].print(tempName);
    }
    // Time Dependent Neumann
    std::cout << "\n--| TIME DEPENDENT NEUMANN:  ";
    std::cout << "\n # TIME DEPENDENT Neumann Bound ID: " << nNeuTimeBound << "\n";
    neuTimeBound.printRowMatlab(" TIME DEPENDENT Neumann Bound ID");
    for (int ibound = 0; ibound < nNeuTimeBound; ibound++)
    {
        std::string tempName = std::to_string(ibound+1) + "th TIME DEPENDENT Neumann function";
        neuTimeMeanFlux[ibound].print(tempName);
    }
}
