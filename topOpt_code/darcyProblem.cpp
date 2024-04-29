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
        alphaIn.initialize((*physics).nNodes_v);
        alpha.length = alphaIn.length;
        alpha.P = alphaIn.P;
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
    STREAM::getRowVector(ParameterFile, line, iss, domains_permeability);

    STREAM::getLines(ParameterFile, line, 3);
    // flag Forcing
    STREAM::getValue(ParameterFile, line, iss, flagForcing);
    if (flagForcing != 0 && flagForcing != 1) throw_line("ERROR: flagForcing different from 0 or 1\n");
    // FORCING
    getline(ParameterFile, line);
    statForcing.resize(dim);
    for (int i = 0; i < dim; i++) STREAM::getValue(ParameterFile, line, iss, statForcing[i]);    
    if (flagForcing == 1)
    {
        timeForcing.resize(dim);
        getline(ParameterFile, line);
        for (int i = 0; i < dim; i++) STREAM::getValue(ParameterFile, line, iss, timeForcing[i]);
    }
    else STREAM::getLines(ParameterFile, line, dim+1);
    
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
    neuMeanP.initialize(nNeuBound);
    if (nNeuBound != 0) 
    {
        getline(ParameterFile, line);
        STREAM::getColVector(ParameterFile, line, iss, neuMeanP, nNeuBound);
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
    neuTimeMeanP.resize(nNeuTimeBound);
    getline(ParameterFile, line);
    for (int ibound = 0; ibound < nNeuTimeBound; ibound++)
    {
        STREAM::getNeuTime(ParameterFile, line, iss, neuTimeBound[ibound], nIdNeuTimeCases[ibound], neuTimeMeanP[ibound]);
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
    std::string NodeFilePath_v  = folderPath + "/Nodes_v.txt";
    std::string ElemFilePath_v  = folderPath + "/Elems_v.txt";
    std::string ElemFilePath  = folderPath + "/Elems.txt";
    std::string BoundElemFilePath_v = folderPath + "/BoundElems_v.txt";
    std::string BoundNodeFilePath_v = folderPath + "/BoundNodes_v.txt";
    std::string BoundElemFilePath = folderPath + "/BoundElems.txt";
    std::string ElemGeoIdsFilePath_v = folderPath + "/ElemGeoIds_v.txt";
    (*physics).bound_info_file_path = BoundElemFilePath;
    (*physics).bound_nodes_v_file_path = BoundNodeFilePath_v;

    //-------OPEN STREAMS-------
    std::ifstream  NodeFile;  NodeFile.open(NodeFilePath, std::ios::in | std::ios::binary);  if (!NodeFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  NodeFile_v;  NodeFile_v.open(NodeFilePath_v, std::ios::in | std::ios::binary);  if (!NodeFile_v.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  ElemFile_v;  ElemFile_v.open(ElemFilePath_v, std::ios::in | std::ios::binary);  if (!ElemFile_v.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  ElemFile;  ElemFile.open(ElemFilePath, std::ios::in | std::ios::binary);  if (!ElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  BoundElemFile_v; BoundElemFile_v.open(BoundElemFilePath_v, std::ios::in | std::ios::binary); if (!BoundElemFile_v.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  BoundNodeFile_v; BoundNodeFile_v.open(BoundNodeFilePath_v, std::ios::in | std::ios::binary); if (!BoundNodeFile_v.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  BoundElemFile; BoundElemFile.open(BoundElemFilePath, std::ios::in | std::ios::binary); if (!BoundElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ifstream  ElemGeoIdsFile_v; ElemGeoIdsFile_v.open(ElemGeoIdsFilePath_v, std::ios::in | std::ios::binary); if (!ElemGeoIdsFile_v.is_open()) throw_line("ERROR, can't open input data file");
    
    //-------------------------
    // READ MESH FROM NODE FILE
    //-------------------------
    std::string line;
    std::istringstream iss;
    STREAM::getValue(NodeFile, line, iss, (*physics).nNodes);
    NodeFile.close();
    STREAM::getValue(NodeFile_v, line, iss, (*physics).nNodes_v);
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    (*physics).nDof = dim * (*physics).nNodes;
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
    
    (*physics).coord_v.initialize(nNodes_v, dim);
    //---
    for(int i = 0; i < nNodes_v; i++)
    {
        prec* p = (*physics).coord_v[i];
        STREAM::getRowVector(NodeFile_v, line, iss, p, dim);
    }
    NodeFile_v.close();

    (*physics).coord.nRow = (*physics).nNodes; (*physics).coord.nCol = dim;
    (*physics).coord.P = (*physics).coord_v.P; (*physics).coord.PP = (*physics).coord_v.PP;
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
    //---
    // velocity elements
    //---
    STREAM::getValue(ElemFile_v, line, iss, (*physics).nElem_v);
    int nElem_v = (*physics).nElem_v;
    (*physics).elem_v.initialize(nElem_v, NodesxElem);
    //---
    for(int i = 0; i < nElem_v; i++)
    {
        int* p = (*physics).elem_v[i];
        STREAM::getRowVector(ElemFile_v, line, iss, p, NodesxElem);
    }
    ElemFile_v.close();

    //--------------------------------
    // READ BOUND INFO FROM BOUND FILE
    //--------------------------------
    std::getline(BoundElemFile_v, line);
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
        std::getline(BoundElemFile_v, line);
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
                STREAM::getRowVector(BoundElemFile_v, line, iss, p, NodesxBoundElem);
                countTime += NodesxBoundElem;
                countElTime++;
            }
            nNeuTimeBoundIdElems[countNeuTime] = currNEl;
            countNeuTime++;
            std::getline(BoundElemFile_v, line);
        } 
        else if (neuBound.hasIn(currID))
        {
            neuBoundIdElems[countNeu] = &neuBoundIdElems_buff[countEl];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &neuBoundIdElems_buff_buff[count];
                neuBoundIdElems_buff[countEl] = p;

                STREAM::getRowVector(BoundElemFile_v, line, iss, p, NodesxBoundElem);
                count += NodesxBoundElem;
                countEl++;
            }
            nNeuBoundIdElems[countNeu] = currNEl;
            countNeu++;
            std::getline(BoundElemFile_v, line);
        } 
        else if (symmBound.hasIn(currID))
        {
            symmBoundIdElems[countSymm] = &symmBoundIdElems_buff[countElSymm];
            for (int j = 0; j < currNEl; j++)
            {
                int* p = &symmBoundIdElems_buff_buff[countS];
                symmBoundIdElems_buff[countElSymm] = p;
                STREAM::getRowVector(BoundElemFile_v, line, iss, p, NodesxBoundElem);
                countS += NodesxBoundElem;
                countElSymm++;
            }
            nSymmBoundIdElems[countSymm] = currNEl;
            countSymm++;
            std::getline(BoundElemFile_v, line);
        } 
        else
        {
            STREAM::getLines(BoundElemFile_v, line, currNEl+1);
            continue;
        }
    }
    // neuBoundIdElems_buff = (int**) realloc(neuBoundIdElems_buff, countEl*sizeof(int*));
    // neuBoundIdElems_buff_buff = (int*) realloc(neuBoundIdElems_buff_buff, count*sizeof(int));
    BoundElemFile_v.close();
    //----- alloc boundIdNodes ---
    std::getline(BoundNodeFile_v, line);
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
        std::getline(BoundNodeFile_v, line);
        iss.str(line);
        iss >> currID; iss >> currNNodes;
        iss.clear();
        boundIdNodes[i] = &boundIdNodes_buff[countNodes];
        for (int j = 0; j < currNNodes; j++)
        {
            STREAM::getValue(BoundNodeFile_v, line, iss, boundIdNodes_buff[countNodes]);
            countNodes++;
        }
        std::getline(BoundNodeFile_v, line);
        nBoundIdNodes[i] = currNNodes;
    }
    BoundNodeFile_v.close();

    //--------------------------------------------
    // RECOVER PRESSURE BOUNDARY ELEMS FOR THE SYMMETRY CONDITIONS
    //--------------------------------------------
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

    int n_geo_entities_ids_v;
    STREAM::getValue(ElemGeoIdsFile_v, line, iss, n_geo_entities_ids_v);
    if (n_geo_entities_ids_v != nElem_v) throw_line("ERROR: nElem_v and n_geo_entities_ids_v are different.\n");
    (*physics).elem_geo_entities_ids_v.initialize(nElem_v);
    STREAM::getColVector(ElemGeoIdsFile_v, line, iss, (*physics).elem_geo_entities_ids_v, nElem_v); 
    ElemGeoIdsFile_v.close();
    (*physics).max_geo_entity_id = (*physics).elem_geo_entities_ids_v.max() + 1;

    //--------------------------------------------
    // ALLOC EVERYTHING USEFUL FOR THE PROBLEM
    //--------------------------------------------
    nDof = dim*nNodes_v + (*physics).nNodes;
    rhs_statForcing.setZeros(nDof); 
    rhs_timeForcing.initialize(nDof); 
    rhs.initialize(nDof);
}

void PROBLEM_DARCY::setBC()
{
    std::cout << "\n Setting Darcy BC\n";
    pause();
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
    std::cout << "\n 1\n";
    pause();
    
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
    std::cout << "\n 2\n";
    pause();

    //--- ITERATE ON THE NO-FLUX BOUNDARIES ----
    int sum = 0; 
    // std::cout << "\n nbound " << nNoFluxBound << "\n";
    // noFluxBound.print();
    // noFluxNod.print();
    // VECTOR_INT::print(nBoundIdNodes, nBound);
    // std::cout << "\n sum " << sum << "\n";
    for (int inoFlux = 0; inoFlux < nNoFluxBound; inoFlux++)
    {
        int ibound = noFluxBound[inoFlux];
        // std::cout << "\n bound: " << ibound << "\n";
        sum += nBoundIdNodes[ibound];
        // std::cout << "\n sum " << sum << "\n";
    }
    std::cout << "\n final sum " << sum << "\n";
    pause();
    noFluxNod.setZeros(sum); //initialize noFluxNod with max possible dim
    std::cout << "\n final sum " << sum << "\n";
    pause();
    nNoFlux = 0;
    
    std::cout << "\n 3\n";
    pause();

    noFluxIdCount.initialize(nNoFluxBound);
    for (int inoFlux = 0; inoFlux < nNoFluxBound; inoFlux++)
    {
        int ibound = noFluxBound[inoFlux];
        for (int inode = 0; inode < nBoundIdNodes[ibound]; inode++)
        {
            int tempNode = boundIdNodes[ibound][inode];
            
            if (boundInfoMat[0][tempNode] == 0)
            {
                boundInfoMat[0][tempNode] = 5;
                noFluxNod[nNoFlux] = tempNode;
                nNoFlux++;
            }
        }
        noFluxIdCount[inoFlux] = nNoFlux;
    }
    noFluxNod.shrink(nNoFlux);
    noFluxNod.length = nNoFlux;
    std::cout << "\n 4\n";
    pause();

    //--- ITERATE ON THE TIME DIR BOUNDARIES ----
    dirTimeIdCount.initialize(nDirTimeBound);
    sum = 0; 
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
    dirTimeVal.initialize(nTimeDir, dim);
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
    dirVal.initialize(nDir, dim);


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
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
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
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[tempElem[inode]][icomp]; 
            }
            prec tempBaseFlux = getArea(matCoord, dim)/dim;
            getNormal(matCoord, normal);

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
                    neuTimeNormal.modifyRow(nTimeNeu, normal);
                    nTimeNeu++;
                }
            };
        }
        neuTimeIdCount[ineu] = nTimeNeu;
    }
    neuTimeNod.shrink(nTimeNeu); neuTimeBaseFlux.shrink(nTimeNeu);
    neuTimeNod.length = nTimeNeu; neuTimeBaseFlux.length = nTimeNeu;
    neuTimeVal.zeros(nTimeNeu, dim); neuTimeNormal.shrinkRows(nTimeNeu);

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
    neuVal.zeros(nNeu, dim);

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

    //--- ITERATE ON THE noFlux BOUNDARIES ----
    sum = 0;
    std::cout << "\n Start Setting\n";
    pause();
    setStatDirBC();
    std::cout << "\n Set Static Dir\n";
    pause();
    setStatNeuBC(boundInfoMat);
    std::cout << "\n Set Static Neu\n";
    pause();
    setSymmBC(boundInfoMat);
    std::cout << "\n Set Symm\n";
    pause();
}
//----------------------------------------------------------------------------
// SET BOUNDARY CONDITIONS
//-----------------------------------------------------------------------------
void PROBLEM_DARCY::setTimeDirBC(prec trialTime) //priority 4
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET TIME DIR BC |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
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
            prec* point  = coord_v[iglob];
            for (int icomp = 0; icomp < dim; icomp++) tempPoint[icomp] = point[icomp];
            tempCount++;
        }

        std::vector<VECTOR> tempTimeDirVal = dirTimeFunc[iDirBound](tempCoord, trialTime);

        tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                dirTimeVal[idir][icomp] = tempTimeDirVal[icomp][tempCount];
            }
            tempCount++;
        }
        old = dirTimeIdCount[iDirBound];
    }
}
//---
void PROBLEM_DARCY::setStatDirBC() //priority 3
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET STAT DIR BC  |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
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
            prec* point  = coord_v[iglob];
            for (int icomp = 0; icomp < dim; icomp++) tempPoint[icomp] = point[icomp];
            tempCount++;
        }

        std::vector<VECTOR> tempDirVal(dim);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            FUNCTION_PARSER::evalCoord(dirFunc[iDirBound][icomp], tempCoord, 0.0, tempDirVal[icomp]);
        }

        tempCount = 0;
        for (int idir = old; idir < endPos; idir++)
        {
            prec* tempDir = dirVal[idir];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                tempDir[icomp] = tempDirVal[icomp][tempCount];
            }

            tempCount++;
        }
        old = endPos;
    }
}
//---
void PROBLEM_DARCY::setTimeNeuBC(prec trialTime) // priority 2
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET TIME NEU BC |--\n----------\n";
    int old = 0;
    prec rho = (*physics).rho;
    int dim = (*physics).dim;
    int neuCount = 0;
    for (int iNeuBound = 0; iNeuBound < nNeuTimeBound; iNeuBound++)
    {
        prec meanP = *neuTimeMeanP[iNeuBound](0, 0, 0, trialTime);
        for (int ineu = old; ineu < neuTimeIdCount[iNeuBound]; ineu++)
        {
            prec* normal = neuTimeNormal[neuCount];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                neuTimeVal[ineu][icomp] = neuTimeBaseFlux[neuCount] * meanP * normal[icomp] / rho;
            }

            neuCount++;
        }
        old = neuTimeIdCount[iNeuBound];
    }
}
//---
void PROBLEM_DARCY::setStatNeuBC(MATRIX_INT& boundInfoMat) // priority 1
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET STAT NEU BC |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
    prec rho = (*physics).rho;
    int* priority = boundInfoMat[0];
    int* neuNodPos = boundInfoMat[1];
    MATRIX matCoord(dim,dim); 
    VECTOR normal(dim);

    for (int ibound = 0; ibound < nNeuBound; ibound++)
    {
        prec meanPByDim = neuMeanP[ibound]/dim;
        //-------------------------------------
        for (int iel = 0; iel < nNeuBoundIdElems[ibound]; iel++)
        {
            int* el = neuBoundIdElems[ibound][iel];
            for (int inode = 0; inode < dim; inode++)
            {
                int iglob = el[inode];
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
            }
            prec Area = getArea(matCoord, dim);
            getNormal(matCoord, normal);
            //---
            prec tempFactor = Area * meanPByDim;
            for (int iloc = 0; iloc < dim; iloc++)
            {
                int iglob = el[iloc];
                // std::cout << "iglob: " << iglob << "\n";
                if (priority[iglob] != 1) 
                {
                    continue;
                }
                int tempNeuNod = neuNodPos[iglob];
                VECTOR value = normal * tempFactor;
                prec* tempNeu = neuVal[tempNeuNod];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    tempNeu[icomp] += value[icomp]/rho;
                }
            }    
        }
    }
}
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
// SET SYMMETRY BC
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
void PROBLEM_DARCY::setSymmBC(MATRIX_INT& boundInfoMat) //priority -1
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS SET SYMMETRY BC |--\n----------\n";
    int dim = (*physics).dim;
    std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;

    int* priority = boundInfoMat[0];
    int* symmNodPos = boundInfoMat[1];

    int nNodes_v = (*physics).nNodes_v;
    int nNodes = (*physics).nNodes;
    prec rho = (*physics).rho;
    prec ni = (*physics).ni;

    MATRIX_INT globElem_v = (*physics).elem_v; int nGlobElem_v = (*physics).nElem_v;
    MATRIX matCoord(dim,dim); 

    std::vector<VECTOR> boundNormal(nSymmBound);
    VECTOR normal(dim);


    VECTOR Volume_v = (*physics).Volume_v;

    std::shared_ptr<prec*[]> Bloc_v = (*physics).Bloc_v.PP;
    std::shared_ptr<prec*[]> Cloc_v = (*physics).Cloc_v.PP;
    std::shared_ptr<prec*[]> Dloc_v = (*physics).Dloc_v.PP;

    CSRMAT JP1; CSRMAT JP2; CSRMAT JP3;
    VECTOR_INT iglobJP1(0); VECTOR_INT kglobJP1(0); VECTOR coefJP1(0);
    VECTOR_INT iglobJP2(0); VECTOR_INT kglobJP2(0); VECTOR coefJP2(0);
    VECTOR_INT iglobJP3(0); VECTOR_INT kglobJP3(0); VECTOR coefJP3(0);

    CSRMAT JV11; CSRMAT JV12; CSRMAT JV13;
    CSRMAT JV21; CSRMAT JV22; CSRMAT JV23;
    CSRMAT JV31; CSRMAT JV32; CSRMAT JV33;
    VECTOR_INT iglobJV11(0); VECTOR_INT jglobJV11(0); VECTOR coefJV11(0);
    VECTOR_INT iglobJV12(0); VECTOR_INT jglobJV12(0); VECTOR coefJV12(0);
    VECTOR_INT iglobJV13(0); VECTOR_INT jglobJV13(0); VECTOR coefJV13(0);
    //
    VECTOR_INT iglobJV21(0); VECTOR_INT jglobJV21(0); VECTOR coefJV21(0);
    VECTOR_INT iglobJV22(0); VECTOR_INT jglobJV22(0); VECTOR coefJV22(0);
    VECTOR_INT iglobJV23(0); VECTOR_INT jglobJV23(0); VECTOR coefJV23(0);
    //
    VECTOR_INT iglobJV31(0); VECTOR_INT jglobJV31(0); VECTOR coefJV31(0);
    VECTOR_INT iglobJV32(0); VECTOR_INT jglobJV32(0); VECTOR coefJV32(0);
    VECTOR_INT iglobJV33(0); VECTOR_INT jglobJV33(0); VECTOR coefJV33(0);


    for (int ibound = 0; ibound < nSymmBound; ibound++)
    {
        //--------------------------------------------------
        int internalElems = (2*dim-2);

        int nElem_v = nSymmBoundIdElems[ibound];

        int nElem = nElem_v/internalElems;

        int** surfElem = symmPressureBoundIdElems[ibound];

        int** surfElem_v = symmBoundIdElems[ibound];

        int nEvalsB = nElem*(dim+1)*internalElems*(dim+1);
        //---------------------------------------------------
        // COLLECT THE ELEMENTS
        //---------------------------------------------------
        int bound = symmBound[ibound];
        VECTOR_INT nodesInBound(nBoundIdNodes[bound]);
        for (int inode = 0; inode < nBoundIdNodes[bound]; inode++)
        {
            nodesInBound[inode] = boundIdNodes[bound][inode];
        }
        //*-*-*-*-*-*-*-*-*-*-*
        int countEl = 0;
        MATRIX_INT newSurfElem_v(nElem_v, dim);
        std::shared_ptr<int[]> globElemFromSurf_v(new int[nElem_v]);
        MATRIX_INT locOrder(nElem_v, dim);
        for (int iel = 0; iel < nGlobElem_v; iel++)
        {
            int* tempElem = globElem_v[iel];
            int count = 0;
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                if (nodesInBound.hasIn(tempElem[iloc]))
                {
                    newSurfElem_v[countEl][count] = tempElem[iloc];
                    locOrder[countEl][count] = iloc;
                    count++;
                } 
            }
            if (count == dim) 
            {
                globElemFromSurf_v[countEl] = iel;
                countEl++;
                if (countEl == nElem_v) break;
            }
        }

        // for (int i = 0; i < nElem_v; i++) 
        // {
        //     for (int icomp = 0; icomp < dim; icomp++) printf("%d ", newSurfElem_v[i][icomp]);

        //     printf("            ");
        //     for (int icomp = 0; icomp < dim; icomp++) printf("%d ", locOrder[i][icomp]);
        //     printf("iel: %d\n", globElemFromSurf_v[i]);
        // }

        //------------------
        // BUILD Jp matrices
        //------------------

        //----------------------
        if (dim == 3)
        {
            
            std::shared_ptr<int[]> i_sparse = VECTOR_INT::makePointer(nEvalsB);
            std::shared_ptr<int[]> k_sparse = VECTOR_INT::makePointer(nEvalsB);
            std::shared_ptr<prec[]> coefJp = VECTOR::makePointer(nEvalsB);

            int count;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                
                count = 0;

                for (int iel = 0; iel < nElem; iel++)
                {
                    int* tempPressureElem = surfElem[iel];
                    // GET NORMAL AND AREA
                    for (int inode = 0; inode < dim; inode++)
                    {
                        int iglob = tempPressureElem[inode];
                        for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
                    }
                    prec Area = getArea(matCoord, dim);
                    getNormal(matCoord, normal);

                    boundNormal[ibound] = normal;

                    prec tempFactor = normal[icomp]*Area/4/rho; // area_v just area/4
                    //---------------------------------------
                    for (int kloc = 0; kloc < 3; kloc++)
                    {
                        int kglob = tempPressureElem[kloc];

                        int miniElementGlobId = 4*iel + kloc;
                        int* miniElement = surfElem_v[miniElementGlobId];

                        for (int iel_v = 0; iel_v < 4; iel_v++)
                        {
                            int globElem_v = 4*iel + iel_v;
                            int* tempVelElem = surfElem_v[globElem_v];

                            for (int iloc = 0; iloc < 3; iloc++)
                            {
                                int iglob = tempVelElem[iloc];
                                
                                prec tempValue = 0;
                                //-----------------------------
                                if (kloc == iel_v)
                                {
                                    if (kglob == iglob) tempValue += 1.0/6 + 1.0/12; // 1/6 + 0.5*2*1/12
                                    else tempValue += 1.0/12 + 0.5*(1.0/12 + 1.0/6);
                                } 
                                //-------------------
                                else if (iel_v == 3)
                                {
                                    int id1 = (kloc+1)%3; int id2 = (kloc+2)%3;
                                    if (iglob == miniElement[id1] || iglob == miniElement[id2]) tempValue = 1.0/8; //0.5*(1/12+1/6)
                                    else tempValue = 1.0/12; // 0.5*(1/12+1/12)
                                } 
                                //-------------------
                                else
                                {
                                    int id1 = (kloc+1)%3; int id2 = (kloc+2)%3;
                                    if (iglob == miniElement[id1] || iglob == miniElement[id2]) tempValue = 1.0/12;
                                    else tempValue = 1.0/24; 
                                }
                                // update
                                i_sparse[count] = iglob;
                                k_sparse[count] = kglob;
                                coefJp[count] = tempValue*tempFactor;
                                count++;
                            }
                        }
                    }
                }  
                switch (icomp)
                {
                    case 0:
                    {
                        iglobJP1.append(i_sparse, count); kglobJP1.append(k_sparse, count); coefJP1.append(coefJp, count);
                        break;
                    }
                    case 1:
                    {
                        iglobJP2.append(i_sparse, count); kglobJP2.append(k_sparse, count); coefJP2.append(coefJp, count);
                        break;
                    }
                    case 2:
                    {
                        iglobJP3.append(i_sparse, count); kglobJP3.append(k_sparse, count); coefJP3.append(coefJp, count);
                        break;
                    }
                }
            }            
            
        } else
        {
            for (int icomp = 0; icomp < dim; icomp++)
            {
                std::shared_ptr<int[]> i_sparse = VECTOR_INT::makePointer(nEvalsB);
                std::shared_ptr<int[]> k_sparse = VECTOR_INT::makePointer(nEvalsB);
                std::shared_ptr<prec[]> coefJp = VECTOR::makePointer(nEvalsB);
                
                int count = 0;
                for (int iel = 0; iel < nElem; iel++)
                {
                    int* tempPressureElem = surfElem[iel];
                    // GET NORMAL AND AREA
                    for (int inode = 0; inode < dim; inode++)
                    {
                        int iglob = tempPressureElem[inode];
                        for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
                    }
                    prec Area = getArea(matCoord, dim);
                    getNormal(matCoord, normal);

                    boundNormal[ibound] = normal;

                    prec tempFactor = normal[icomp]*Area/2/rho; // length_v just length/2
                    //---------------------------------------
                    for (int kloc = 0; kloc < 2; kloc++)
                    {
                        int kglob = tempPressureElem[kloc];

                        for (int iel_v = 0; iel_v < 2; iel_v++)
                        {
                            int globElem_v = 2*iel + iel_v;
                            int* tempVelElem = surfElem_v[globElem_v];

                            for (int iloc = 0; iloc < 2; iloc++)
                            {
                                int iglob = tempVelElem[iloc];
                                prec tempValue = 0;
                                if (iloc == kloc) tempValue = 1.0/3+1.0/12; // 1/3 + 0.5*1/6
                                else tempValue = 1.0/3; // 1/6+0.5*1/3
                                // update
                                i_sparse[count] = iglob;
                                k_sparse[count] = kglob;
                                coefJp[count] = tempValue*tempFactor;
                                count++;
                            }
                        }
                    }
                }
                // BUILD LOCAL CSRMAT 
                // ADD CSRMAT TO VECTOR
                switch (icomp)
                {
                    case 0:
                    {
                        iglobJP1.append(i_sparse, count); kglobJP1.append(k_sparse, count); coefJP1.append(coefJp, count);
                        break;
                    }
                    case 1:
                    {
                        iglobJP2.append(i_sparse, count); kglobJP2.append(k_sparse, count); coefJP2.append(coefJp, count);
                        break;
                    }
                }
            }
        }
        //------------------
        // BUILD Jv matrices
        //------------------
        int nEvals = nElem_v*dim*dim;

        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                std::shared_ptr<int[]> i_sparse = VECTOR_INT::makePointer(nEvals);
                std::shared_ptr<int[]> j_sparse = VECTOR_INT::makePointer(nEvals);
                std::shared_ptr<prec[]> coefJv = VECTOR::makePointer(nEvals);
                int count = 0;
                
                for (int iel = 0; iel < nElem_v; iel++)
                {
                    int* tempElem = newSurfElem_v[iel];
                    // GET NORMAL AND AREA
                    for (int inode = 0; inode < dim; inode++)
                    {
                        int iglob = tempElem[inode];
                        for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
                    }
                    prec Area = getArea(matCoord, dim);
                    

                //    getNormal(matCoord, normal);

                normal = boundNormal[ibound];

                    
                    int globEl = globElemFromSurf_v[iel];

                    prec tempFactor = -ni*Area/dim*normal[icomp]*normal[jcomp];  //THE MINUS IS HERE!!!

                    // if (icomp != jcomp) tempFactor = 0;
                    for (int iloc = 0; iloc < dim; iloc++)
                    {
                        int iglob = tempElem[iloc];

                        for (int jloc = 0; jloc < dim; jloc++)
                        {
                            int jglob = tempElem[jloc];

                            int locIndex = locOrder[iel][jloc];
                            prec b_j = Bloc_v[globEl][locIndex]; prec c_j = Cloc_v[globEl][locIndex]; 
                            prec tempValue = normal[0]*b_j + normal[1]*c_j;
                            if (dim == 3)
                            {
                                prec d_j = Dloc_v[globEl][locIndex]; 
                                tempValue += normal[2]*d_j;
                            } 
                            //---
                            i_sparse[count] = iglob;
                            j_sparse[count] = jglob; 
                            coefJv[count]   = tempValue*tempFactor;
                            count++;
                            //---
                        }
                    }
                }
                switch (icomp)
                {
                    case 0:
                    {
                        switch (jcomp)
                        {
                            case 0:
                            {
                                iglobJV11.append(i_sparse, count); jglobJV11.append(j_sparse, count); coefJV11.append(coefJv, count);
                                break;
                            }
                            case 1:
                            {
                                iglobJV12.append(i_sparse, count); jglobJV12.append(j_sparse, count); coefJV12.append(coefJv, count);
                                break;
                            }
                            case 2:
                            {
                                iglobJV13.append(i_sparse, count); jglobJV13.append(j_sparse, count); coefJV13.append(coefJv, count);
                                break;
                            }
                        }
                        
                        break;
                    }
                    case 1:
                    {
                        switch (jcomp)
                        {
                            case 0:
                            {
                                iglobJV21.append(i_sparse, count); jglobJV21.append(j_sparse, count); coefJV21.append(coefJv, count);
                                break;
                            }
                            case 1:
                            {
                                iglobJV22.append(i_sparse, count); jglobJV22.append(j_sparse, count); coefJV22.append(coefJv, count);
                                break;
                            }
                            case 2:
                            {
                                iglobJV23.append(i_sparse, count); jglobJV23.append(j_sparse, count); coefJV23.append(coefJv, count);
                                break;
                            }
                        }
                        break;
                    }
                    case 2:
                    {
                        switch (jcomp)
                        {
                            case 0:
                            {
                                iglobJV31.append(i_sparse, count); jglobJV31.append(j_sparse, count); coefJV31.append(coefJv, count);
                                break;
                            }
                            case 1:
                            {
                                iglobJV32.append(i_sparse, count); jglobJV32.append(j_sparse, count); coefJV32.append(coefJv, count);
                                break;
                            }
                            case 2:
                            {
                                iglobJV33.append(i_sparse, count); jglobJV33.append(j_sparse, count); coefJV33.append(coefJv, count);
                                break;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
    if (nSymmBound > 0)
    {
        //------------------------------
        // BUILD LOCAL CSRMAT JP
        int globCount = iglobJP1.length;
        JP1.initialize(nNodes_v, nNodes, globCount, iglobJP1.P, kglobJP1.P, coefJP1.P);
        globCount = iglobJP2.length;
        JP2.initialize(nNodes_v, nNodes, globCount, iglobJP2.P, kglobJP2.P, coefJP2.P);

        if (dim == 3)
        {
            globCount = iglobJP3.length;
            JP3.initialize(nNodes_v, nNodes, globCount, iglobJP3.P, kglobJP3.P, coefJP3.P);
        }
        //------------------------
        // BUILD LOCAL CSRMAT JV
        globCount = iglobJV11.length;
        JV11.initialize(nNodes_v, nNodes_v, globCount, iglobJV11.P, jglobJV11.P, coefJV11.P);
        globCount = iglobJV12.length;
        JV12.initialize(nNodes_v, nNodes_v, globCount, iglobJV12.P, jglobJV12.P, coefJV12.P);
        //
        globCount = iglobJV21.length;
        JV21.initialize(nNodes_v, nNodes_v, globCount, iglobJV21.P, jglobJV21.P, coefJV21.P);
        globCount = iglobJV22.length;
        JV22.initialize(nNodes_v, nNodes_v, globCount, iglobJV22.P, jglobJV22.P, coefJV22.P);
        //-----------------------
        std::vector<std::vector<CSRMAT*>> tempVec(dim);

        for (int icomp = 0; icomp < dim; icomp++)
        {
            tempVec[icomp].resize(dim+1);
        }
        //-----------------------
        if (dim == 3)
        {
            globCount = iglobJP3.length;
            JP3.initialize(nNodes_v, nNodes, globCount, iglobJP3.P, kglobJP3.P, coefJP3.P);


            globCount = iglobJV13.length;
            JV13.initialize(nNodes_v, nNodes_v, globCount, iglobJV13.P, jglobJV13.P, coefJV13.P);
            globCount = iglobJV23.length;
            JV23.initialize(nNodes_v, nNodes_v, globCount, iglobJV23.P, jglobJV23.P, coefJV23.P);
            //
            globCount = iglobJV31.length;
            JV31.initialize(nNodes_v, nNodes_v, globCount, iglobJV31.P, jglobJV31.P, coefJV31.P);
            globCount = iglobJV32.length;
            JV32.initialize(nNodes_v, nNodes_v, globCount, iglobJV32.P, jglobJV32.P, coefJV32.P);

            globCount = iglobJV33.length;
            JV33.initialize(nNodes_v, nNodes_v, globCount, iglobJV33.P, jglobJV33.P, coefJV33.P);

            //-------------------------------------------
            tempVec[0][0] = &JV11; tempVec[0][1] = &JV12; tempVec[0][2] = &JV13; tempVec[0][3] = &JP1;
            tempVec[1][0] = &JV21; tempVec[1][1] = &JV22; tempVec[1][2] = &JV23; tempVec[1][3] = &JP2;
            tempVec[2][0] = &JV31; tempVec[2][1] = &JV12; tempVec[2][2] = &JV33; tempVec[2][3] = &JP3;
        }
        else
        {
            tempVec[0][0] = &JV11; tempVec[0][1] = &JV12; tempVec[0][2] = &JP1;
            tempVec[1][0] = &JV21; tempVec[1][1] = &JV22; tempVec[1][2] = &JP2;            
        }

        // int L = nNodes_v*dim*dim; 
        // std::shared_ptr<int[]> i_sparse(new int[L]);
        // std::shared_ptr<int[]> j_sparse(new int[L]);
        // std::shared_ptr<prec[]> coef_sparse(new prec[L]);
        // int countTemp = 0;
        // for (int icomp = 0; icomp < dim; icomp++)
        // {
        //     for (int jcomp = 0; jcomp < dim; jcomp++)
        //     {
        //         for (int i = 0; i < nNodes_v; i++)
        //         {
        //             i_sparse[countTemp] = i + nNodes_v*icomp;
        //             j_sparse[countTemp] = i + nNodes_v*jcomp;
        //             coef_sparse[countTemp] = 0;
        //             countTemp++;
        //         }
        //     }
        // }

        // J.initialize(nDof, nDof, L, i_sparse, j_sparse, coef_sparse);

        PARALLEL::createBlockMatrix(tempVec, J);

        J.enlargeRows(nDof);

        std::shared_ptr<int[]> iat = J.iat;
        std::shared_ptr<prec[]> coef = J.coef;
        for (int ineu = 0; ineu < nNeu; ineu++)
        {
            int iglob = neuNod[ineu];

            for (int icomp = 0; icomp < dim; icomp++)
            {
                int globRow = icomp*nNodes_v + iglob;
                int startPos = iat[globRow]; int endPos = iat[globRow+1];

                for (int pos = startPos; pos < endPos; pos++)
                {
                    coef[pos] = 0;
                }
            }
        }
        for (int ineu = 0; ineu < nTimeNeu; ineu++)
        {
            int iglob = neuTimeNod[ineu];

            for (int icomp = 0; icomp < dim; icomp++)
            {
                int globRow = icomp*nNodes_v + iglob;
                int startPos = iat[globRow]; int endPos = iat[globRow+1];

                for (int pos = startPos; pos < endPos; pos++)
                {
                    coef[pos] = 0;
                }
            }
        }
    }
    // PREPARE FOR U*n = 0

    VECTOR_INT tempCounter(nSymm); PARALLEL::resetZeros(tempCounter);
    VECTOR_INT passed(nSymm); passed.reset(-1.0);

    int countDiscard = 0;
    for (int ibound = 0; ibound < nSymmBound; ibound++)
    {
        int countEstup = 0;

        int nElem_v = nSymmBoundIdElems[ibound];

        int** surfElem_v = symmBoundIdElems[ibound];

        VECTOR_INT tempSymmPos(nNodes_v); int* tempSymmPosP = &(tempSymmPos[0]); 
        PARALLEL::copy(symmNodPos, nNodes_v, tempSymmPosP);

        for (int iel = 0; iel < nElem_v; iel++)
        {
            int* tempElem = surfElem_v[iel];
            // GET NORMAL
            for (int inode = 0; inode < dim; inode++)
            {
                int iglob = tempElem[inode];
                for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord_v[iglob][icomp]; 
            }
            getNormal(matCoord, normal);

            // normal = boundNormal[ibound];

            for (int iloc = 0; iloc < dim; iloc++)
            {
                int iglob = tempElem[iloc];
                if (priority[iglob] != -1) 
                {
                    continue;
                }
                //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                int isymm = tempSymmPos[iglob];
                
                if (passed[isymm] == -1)
                {
                    passed[isymm] = ibound; 
                }
                //
                if (passed[isymm] != ibound)
                {
                    tempSymmPos[iglob] = nSymm;
                    isymm = nSymm;
                    nSymm++;
                    passed.enlarge(nSymm); passed[nSymm-1] = ibound;
                    tempCounter.enlarge(nSymm); tempCounter[nSymm-1] = 0;
                    symmVal.enlargeRows(nSymm);

                    symmNod.enlarge(nSymm); symmNod[nSymm-1] = iglob;
                    symmRealIglob.enlarge(nSymm);

                    countDiscard++;

                    countEstup++;
                }

                prec* tempSymm = symmVal[isymm];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    tempSymm[icomp] = normal[icomp];
                }
                tempCounter[isymm]++;
            }
        }
    }

    // for (int iSymm = 0; iSymm < nSymm; iSymm++)
    // {
    //     int tempCount = tempCounter[iSymm];
    //     for (int icomp = 0; icomp < dim; icomp++)
    //     {
    //         symmVal[iSymm][icomp] /= tempCount;
    //     }
    // }
}
//-----------------------------------------------------------------------
// FORCING 
//-----------------------------------------------------------------------
void PROBLEM_DARCY::applyStatForcing()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS APPLY STATIC FORCING |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    int nElem_v = (*physics).nElem_v;

    std::vector<VECTOR> forcingAtNodes(dim);
    for (int icomp = 0; icomp < dim; icomp++)
    {
        forcingAtNodes[icomp].initialize(nNodes_v);
        FUNCTION_PARSER::evalCoord(statForcing[icomp], (*physics).coord_v, 0, forcingAtNodes[icomp]);
    }
    //----------
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int* tempElem = elem_v[iel];
        for (int icomp = 0; icomp < dim; icomp++)
        {
            prec tempForcingSum = 0;
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = tempElem[iloc];
                tempForcingSum += forcingAtNodes[icomp][iglob];
            }
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = tempElem[iloc];
                prec forcingFactor = (forcingAtNodes[icomp][iglob] + tempForcingSum)/(dim+2);
                int globId = icomp*nNodes_v+iglob;
                rhs_statForcing[globId] += Volume_v[iel]/(dim+1)*forcingFactor;
            }
        }
    }        
    //---------------------------------------
}
//---
void PROBLEM_DARCY::updateTimeForcing()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS UPDATE TIME FORCING |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;
    int nElem_v = (*physics).nElem_v;
    std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    std::shared_ptr<prec[]> Volume_v = (*physics).Volume_v.P;
    prec evalTime = time + deltaT;
    rhs_timeForcing.resetZeros();
    std::vector<VECTOR> forcingAtNodes(dim);

    for (int icomp = 0; icomp < dim; icomp++)
    {    
        
        forcingAtNodes[icomp].initialize(nNodes_v);
        FUNCTION_PARSER::evalCoord(timeForcing[icomp], (*physics).coord_v, evalTime, forcingAtNodes[icomp]);
        
        //std::cout << "\n" << forcingAtNodes[icomp].length << "\n";
    }
    //----------
    
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int* tempElem = elem_v[iel];
        
        for (int icomp = 0; icomp < dim; icomp++)
        {
            prec tempForcingSum = 0;
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                
                int iglob = tempElem[iloc];                
                tempForcingSum += forcingAtNodes[icomp][iglob];
                
            }
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = tempElem[iloc];
                prec forcingFactor = (forcingAtNodes[icomp][iglob] + tempForcingSum)/(dim+2);
                int globId = icomp*nNodes_v+iglob;
                rhs_timeForcing[globId] += Volume_v[iel]/(dim+1)*forcingFactor;
            }
            
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
    imposeTimeBC(SYSMAT_final, rhs);
}

//------------------------------
// IMPOSE STATIC BC
//------------------------------
void PROBLEM_DARCY::imposeStaticBC(CSRMAT &SYSMAT_NS, VECTOR &rhs)
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
            // noFlux dir
            for (int inoFlux = 0; inoFlux < nNoFlux; inoFlux++)
            {
                int iglob = noFluxNod[inoFlux] + nNodes_v*icomp; 
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
            // noFlux
            // set to zero row & col entries
            int countPosToZero = 0; int countIglobToZero = 0; int countPosToOne = 0;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int inoFlux = 0; inoFlux < nNoFlux; inoFlux++)
                {
                    int iglob = noFluxNod[inoFlux] + icomp*nNodes_v; 
                    VECTOR_INT rowInColToZero = SYSMAT_NS.setColInRowToZero(iglob, noFluxPosToZero, countPosToZero, noFluxPosToOne, countPosToOne);  
                    rhs[iglob] = 0;                                             
                    for (int i = 0; i < rowInColToZero.length; i++) 
                    {
                        int tempPos = SYSMAT_NS(rowInColToZero[i], iglob, 0);
                        noFluxPosToZero[countPosToZero] = tempPos; countPosToZero++;
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
            //--- noFlux
            int nToZero = noFluxPosToZero.length; int nToOne = noFluxPosToOne.length;
            std::shared_ptr<prec[]> coefMat = SYSMAT_NS.coef;
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int inoFlux = 0; inoFlux < nNoFlux; inoFlux++)
                {
                    int iglob = noFluxNod[inoFlux] + icomp*nNodes_v;
                    rhs[iglob] = 0.0;
                }
            }
            //---
            for (int countPosToZero = 0; countPosToZero < nToZero; countPosToZero++)
            {
                int pos = noFluxPosToZero[countPosToZero];
                coefMat[pos] = 0.0;
            }
            for (int countPosToOne = 0; countPosToOne < nToOne; countPosToOne++)
            {
                int pos = noFluxPosToOne[countPosToOne];
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
void PROBLEM_DARCY::imposeTimeBC(CSRMAT &SYSMAT_NS, VECTOR &rhs)
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
// UPDATE BC
//-------------------------------------------------
void PROBLEM_DARCY::updateBC(prec trialTime)
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS UPDATE BC  |--\n----------\n";
    setTimeDirBC(trialTime);
    setTimeNeuBC(trialTime);
}

//-------------------------------------------------
// UPDATE SYSMAT_NS 
//-------------------------------------------------
void PROBLEM_DARCY::updateSYSMAT()
{
    // if ((*physics).completeLog == 0) std::cout << "\n----------\n--| NS UPDATE SYSMAT_NS |--\n----------\n";
    
    // SYSMAT = SYSMAT_base; // initialize SYSMAT as SYSMAT_base in order to add updated matrices in correct blocks

    // // add Time Derivative Matrix
    // std::shared_ptr<prec[]> coefM = M.coef;
    // prec factorP = 1 / deltaT;
    // // add P
    // addToSysmat(coefM, factorP); // add in diagonal blocks

    // // convective term assembled in the convective convergence loop

    // // add Optimization Variable Matrix
    // // assemble
    // assembleMa();
    // std::shared_ptr<prec[]> coefMa = Ma.P;
    // prec factorMa = 1 / (*physics).rho;
    // // add Ma
    // addToSysmat(coefMa, factorMa); // add in diagonal blocks

}
//---

//-------------------------------------------------
// UPDATE RHS
//-------------------------------------------------
void PROBLEM_DARCY::updateRHS()
{
    if ((*physics).completeLog == 0) std::cout << "\n----------\n--| UPDATE RHS |--\n----------\n";
    int dim = (*physics).dim;
    int nNodes_v = (*physics).nNodes_v;

    if (flagForcing == 1)
    {
        updateTimeForcing();
        VECTOR::sum(rhs_statForcing, rhs_timeForcing, rhs);
        currForcingAtNodes = rhs;
        for (int icomp = 0; icomp < dim; icomp++)      //(P * lastSol)/deltaT;
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
    }
    else
    {
        static bool alreadyVisited = false;
        if (!alreadyVisited)
        {
            alreadyVisited = true;
            currForcingAtNodes = rhs_statForcing;
        } 
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

            prec* rhs_fP = &rhs_statForcing[startId];
            VECTOR::sum(rhs_fP, tSol.P, nNodes_v, tempP);
        }
    }
    
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
    // assemble();
    // applyStatForcing();
    // // prepare BC update in case of lifting functions imposition
    // int dim = (*physics).dim;
    // int nElem_v = (*physics).nElem_v;
    // std::shared_ptr<int*[]> elem_v = (*physics).elem_v.PP;
    // if ((*physics).flagBC == 1)
    // {
    //     std::shared_ptr<int[]> iat = SYSMAT_base.iat;
    //     int countPosToZero = 0;
    //     //--- noFlux
    //     for (int inoFlux = 0; inoFlux < nNoFlux; inoFlux++)
    //     {
    //         int iglob = noFluxNod[inoFlux];
    //         int nnzInRow = iat[iglob+1] - iat[iglob];
    //         countPosToZero += nnzInRow;
    //     }
    //     countPosToZero -= nNoFlux;
    //     countPosToZero *= 2;
    //     noFluxPosToZero.initialize(countPosToZero*dim);
    //     noFluxPosToOne.initialize(nNoFlux*dim);
    //     //--- static
    //     countPosToZero = 0;
    //     for (int idir = 0; idir < nDir; idir++)
    //     {
    //         int iglob = dirNod[idir];
    //         int nnzInRow = iat[iglob+1] - iat[iglob];
    //         countPosToZero += nnzInRow;
    //     }
    //     countPosToZero -= nDir;
    //     dirIglobToZero.initialize(countPosToZero * dim);
    //     countPosToZero *= 2;
    //     dirPosToZero.initialize(countPosToZero * dim);
    //     dirPosToOne.initialize(nDir * dim);
    //     dirCountToZero.initialize(nDir*dim + 1);
    //     //--- time-dependent
    //     countPosToZero = 0;
    //     for (int idir = 0; idir < nTimeDir; idir++)
    //     {
    //         int iglob = dirTimeNod[idir];
    //         int nnzInRow = iat[iglob+1] - iat[iglob];
    //         countPosToZero += nnzInRow;
    //     }
    //     countPosToZero -= nTimeDir;
    //     dirTimeIglobToZero.initialize(countPosToZero * dim);
    //     countPosToZero *= 2;
    //     dirTimePosToZero.initialize(countPosToZero * dim);
    //     dirTimePosToOne.initialize(nTimeDir * dim);
    //     dirTimeCountToZero.initialize(nTimeDir*dim + 1);
    //     //--------------------------------------------------
    //     // SYMMETRY BC
    //     //--------------------------------------------------
    //     countPosToZero = 0;
    //     for (int isymm = 0; isymm < nSymm; isymm++)
    //     {
    //         int iglob = symmNod[isymm];
    //         int nnzInRow = iat[iglob+1] - iat[iglob];
    //         countPosToZero += nnzInRow;
    //     }
    //     countPosToZero -= nSymm*dim;
    //     symmPosToZero.initialize(countPosToZero);
    //     symmPosToCoef.initialize(nSymm*dim);
    // }
    
    // //-----------------------
    // // FREE USELESS STUFF
    // //-----------------------
    // boundIdNodes_buff.reset(); boundIdNodes.reset(); 
    // // dirBound.dlt();
    // // noFluxBound.dlt(); 
    // //----------------------------------------------------------
    // // get max length of each element
    // //----------------------------------------------------------
    // (*physics).h_v.initialize(nElem_v);
    // std::shared_ptr<prec*[]> coord_v = (*physics).coord_v.PP;
    // std::vector<VECTOR> tempCoord(dim+1);
    // for (int iel = 0; iel < nElem_v; iel++)
    // {
    //     prec maxDist = 0;
    //     int* tempElem = elem_v[iel];
    //     //----------------------------------------
    //     for (int iloc = 0; iloc < dim+1; iloc++)
    //     {
    //         int iglob = tempElem[iloc];
    //         prec* nodeCoord = coord_v[iglob];
    //         tempCoord[iloc].initialize(dim);
    //         for (int icomp = 0; icomp < dim; icomp++)
    //         {
    //             tempCoord[iloc][icomp] = nodeCoord[icomp];
    //         }
    //     }
    //     //----------------------------------------
    //     for (int iloc = 0; iloc < dim+1; iloc++)
    //     {
    //         for (int jloc = iloc+1; jloc < dim+1; jloc++)
    //         {
    //             VECTOR diff = tempCoord[iloc] - tempCoord[jloc];
    //             prec tempDist = diff.norm();
    //             if (tempDist > maxDist)
    //             {
    //                 maxDist = tempDist;
    //             }
    //         }
    //     }
    //     (*physics).h_v[iel] = maxDist;
    // }
    // //--------------------------------------
    // // PRINT INIT COND
    // //--------------------------------------
    // int nNodes = (*physics).nNodes;
    // int nElem = (*physics).nElem;
    // int nNodes_v = (*physics).nNodes_v;
    // std::string folderName;
    // switch ((*physics).isNS)
    // {
    //     case 0:
    //     {
    //         folderName = "S_sol";
    //         break;
    //     }
    //     case 1:
    //     {
    //         folderName = "NS_sol";
    //         break;
    //     }
    //     default:
    //     {
    //         break;
    //     }
    // }
    // VTKWriter.initializeForNS(name, dim, nNodes, nElem, folderName);

    // setInitCond();

    // // GET INIT VELOCITY
    // MATRIX coord(nNodes, dim);
    // coord = (*physics).coord;
    // MATRIX_INT elem(nElem,dim+1);
    // elem = (*physics).elem;
    // MATRIX velocity(nNodes, dim);
    // for (int icomp = 0; icomp < dim; icomp++)
    // {
    //     for (int inode = 0; inode < nNodes; inode++)
    //     {
    //         int iglob = inode + nNodes_v*icomp;
    //         velocity[inode][icomp] = lastSol[iglob];
    //     }
    // }
    // // GET INIT PRESSURE
    // VECTOR pressure(nNodes);
    // prec* tempP = &lastSol[dim*nNodes_v];
    // pressure = tempP;

    // // VTKWriter.write(coord, elem, pressure, 0, deltaT, velocity);
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
    throw_line("ERROR: non updated solution case. Disused method\n");
    time = 0;
    // prec t_end = (*physics).t_end;
    deltaT = (*physics).deltaT;
    globIter = 0;

    oneStepSolver(1e-3, 5);
    simulation_times_solution.append_row(lastSol);
    print_sol_in_VTK((*physics).NS_solution);

    VTKWriter.closeTFile();
}
//----

void PROBLEM_DARCY::StatSolverIterative()
{
    
    prec t_end = (*physics).t_end; // = (*physics).solution_times[0]
    deltaT = (*physics).deltaT; // fixed value defined in the geometry.h time settings method
    time = 0;
    globIter = 0;
    while (time < t_end)
    {
        oneStepSolver();
        if (time+deltaT > t_end) 
        {
            deltaT = t_end-time;
            // convergence_toll = 1e-4;
            // convergence_it = 20;
        }
    }
    // std::cout << "\nNS_las_sol_norm: " << lastSol.norm() << "\n";
    // pause();
    (*physics).NS_solution.modifyRow(0, lastSol);
    print_sol_in_VTK((*physics).NS_solution);
    VTKWriter.closeTFile();
}

void PROBLEM_DARCY::Solver()
{
    real_solution_times.resetZeros();
    real_solution_times.initialize(0);
    simulation_times_solution.complete_reset();
    // (*physics).NS_solution.complete_reset();
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
    simulation_times_solution.append_row(init_cond);
    //print_one_step_sol_in_VTK((*physics).dim, (*physics).nNodes, (*physics).nNodes_v, 0);

    while ((t_end-time) > (1e-15 * (nTimeSteps+1)))
    {
        oneStepSolver();
        simulation_times_solution.append_row(lastSol);
        // std::cout << "\nNS_las_sol_norm: " << lastSol.norm() << "\n";
        if (time+deltaT > t_end) deltaT = t_end-time;
    }
    // pause();
    print_sol_in_VTK((*physics).NS_solution);
    VTKWriter.closeTFile();
}

//-------------------------
// ONE-STEP SOLVER
//-------------------------
// General
void PROBLEM_DARCY::oneStepSolver(prec toll, int itMax)
{
    std::cout << "\nCREATE ONE STEP SOLVER\n";
    pause();
}

//-------------------------------
// PRINT ONE STEP SOLUTION IN VTK
//-------------------------------
void PROBLEM_DARCY::print_one_step_sol_in_VTK(int dim, int nNodes, int nNodes_v, prec trial_time)
{
    // // GET VELOCITY
    // MATRIX velocity = getVelocityFromSol(lastSol);
    // // GET PRESSURE
    // VECTOR pressure(nNodes);
    // prec* tempP = &lastSol[dim*nNodes_v];
    // pressure = tempP;
    
    // //-----------------
    // if (printRes) VTKWriter.write((*physics).coord, (*physics).elem, trial_time, pressure, 1, "Pressure", velocity, dim, "Velocity");
}

void PROBLEM_DARCY::print_sol_in_VTK(MATRIX &requested_sol)
{
    // for (int itime = 0; itime < (*physics).solution_times.length; itime++)
    // {
    //     prec trial_time;
    //     if ((*physics).isStationary == 1) //stationary case: only one time = "infinity", so it has no meaning to report it.
    //     {
    //         trial_time = 0;
    //     }
    //     else // time dependent case: starting from t0.
    //     {
    //         trial_time = (*physics).solution_times[itime];
    //     }
         
    //     // GET VELOCITY
    //     VECTOR temp_sol = requested_sol.get_row(itime);
    //     // temp_sol.printRowMatlab("prova");
        
    //     MATRIX velocity = getVelocityFromSol(temp_sol);
    //     // MATRIX::printForMatlab(velocity, "vel");
    //     // std::cout << "tempsolL: " << temp_sol.length << "\n";
    //     // std::cout << "nDof: " << nDof << "\n";
    //     // pause();
    //     // GET PRESSURE
    //     VECTOR pressure((*physics).nNodes);
    //     prec* tempP = &temp_sol[(*physics).dim*(*physics).nNodes_v];
    //     pressure = tempP;
        
    //     //-----------------
    //     if (printRes) VTKWriter.write((*physics).coord, (*physics).elem, trial_time, pressure, 1, "Pressure", velocity, (*physics).dim, "Velocity");
    // }
         
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
    neuMeanP.printRowMatlab(" STATIC Neumann flux");

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
        neuTimeMeanP[ibound].print(tempName);
    }
}
