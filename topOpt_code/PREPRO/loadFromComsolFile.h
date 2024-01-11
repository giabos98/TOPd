#include "../CODE_HEADERS/codeHeader.h"

namespace fs = std::filesystem;

int loadMesh(char* fileName, int N, const char* infolderPath, bool BINREAD)
{
    char c = 'a';
    int first_id = 0;
    int last_id = N;
    while (c != '.' && last_id != 0)
    {
        last_id--;
        c = fileName[last_id]; 
    }
    if (last_id == 0) throw_line("ERROR in load file, missing type specifier, like \".dat\" ");

    const char* prefix = "PROBLEM_DATA/";
    int N_prefix = strlen(prefix);
    char* folderName = (char*) malloc((last_id+N_prefix) * sizeof(char));

    std::string name;
    for (int i = 0; i < N_prefix; i++)
    {
        folderName[i] = prefix[i];
    }

    for (int i = 0; i < last_id; i++)
    {
        int j = i + N_prefix;
        int k = i + first_id;
        folderName[j] = fileName[k];
        name += fileName[k];
    }
    // CREATE NEW FOLDER IN THE "data" FOLDER FOR THE SPECIFIC PROBLEM
    fs::create_directories(folderName);
    //-----------------------------------------------------------------
    // READ COMSOL FILE AND SAVE DATA IN APPROPRIATE ".txt" FILES
    std::ifstream inFile;
    std::string inFilePath = infolderPath;
    inFilePath = inFilePath + fileName;
    if (BINREAD)
    {
        inFile.open(inFilePath, std::ios::binary);
    }
    else inFile.open(inFilePath);
    std::string folderPath = folderName;
    std::string NodeFilePath = folderPath + "/Nodes.txt";
    std::string NodeComSolFilePath = folderPath + "/NodesComSol.txt";
    std::string ElemFilePath = folderPath + "/Elems.txt";
    std::string BoundElemFilePath = folderPath + "/BoundElems.txt";
    std::string BoundNodeFilePath = folderPath + "/BoundNodes.txt";
    std::string ElemGeoIdsFilePath = folderPath + "/ElemGeoIds.txt";
    
    std::ofstream  NodeFile;  NodeFile.open(NodeFilePath, std::ios::out | std::ios::binary);  if (!NodeFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ofstream  NodeComSolFile;  NodeComSolFile.open(NodeComSolFilePath, std::ios::out | std::ios::binary);  if (!NodeComSolFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ofstream  ElemFile;  ElemFile.open(ElemFilePath, std::ios::out | std::ios::binary);  if (!ElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ofstream  BoundElemFile; BoundElemFile.open(BoundElemFilePath, std::ios::out | std::ios::binary); if (!BoundElemFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ofstream  BoundNodeFile; BoundNodeFile.open(BoundNodeFilePath, std::ios::out | std::ios::binary); if (!BoundNodeFile.is_open()) throw_line("ERROR, can't open input data file");
    std::ofstream  ElemGeoIdsFile; ElemGeoIdsFile.open(ElemGeoIdsFilePath, std::ios::out | std::ios::binary); if (!ElemGeoIdsFile.is_open()) throw_line("ERROR, can't open input data file");

    if (!inFile.is_open()) throw_line("ERROR, can't open input data file");
    //----------------------------
    std::string line;
    for (int i = 0; i < 16; i++) std::getline(inFile, line);
    std::getline(inFile, line);
    std::istringstream iss(line);
    //---------------------------
    int nDim;
    iss >> nDim;
    //--------------------------
    // LOAD NODES
    //-------------------------
    std::getline(inFile, line);
    iss.clear();
    iss.str(line);
    int Nodes; iss >> Nodes;
    std::getline(inFile, line); std::getline(inFile, line); std::getline(inFile, line);
    // while(std::getline(inFile, line, '#'))

    
    //printf("-------------------------\nStarting Nodes Load \n-------------------------\n");
    NodeFile << Nodes << " \n";
    // int countPerc = 1; int maxCountPerc = Nodes/100;
    MATRIX coord(Nodes, nDim);
    for(int i = 0; i < Nodes; i++)
    {
        std::getline(inFile, line);
        iss.clear();
        iss.str(line);

        prec* tempNode = coord[i];
        for (int d = 0; d < nDim; d++)
        {
            prec tempCoord_d;
            iss >> tempCoord_d;
            NodeFile << tempCoord_d << " ";
            NodeComSolFile << tempCoord_d << " ";
            tempNode[d] = tempCoord_d;
        }
        NodeFile << "\n";
        NodeComSolFile << "\n";

        //print percentage of loading
        // if (countPerc == maxCountPerc)
        // {
        //     prec tempPerc = (prec) 100*(i+1)/Nodes;
        //     //printf("Loading Nodes: %4.2" format "%% \n", tempPerc);
        //     countPerc = 0;
        // }
        // countPerc++;
    }
    NodeFile.close();
    NodeComSolFile.close();
    //printf("-------------------------\nNodes loaded correctly \n-------------------------\n");

    // skip appropriate number of lines
    for (int i = 0; i < 8; i++) std::getline(inFile, line);
    //--------------------------
    // LOAD BOUNDS
    //-------------------------
    std::getline(inFile, line);
    iss.clear();
    iss.str(line);
    //---------------------------
    int boundNodesxElem;
    int NBoundelem;
    iss >> boundNodesxElem;
    std::getline(inFile, line);
    iss.clear();
    iss.str(line); iss >> NBoundelem;
    std::getline(inFile, line);


    MATRIX_INT BoundElems(NBoundelem, boundNodesxElem);
    //printf("-------------------------\nStarting Boundary Element Read \n-------------------------\n");
    STREAM::getMatrix(inFile, line, iss, BoundElems);
    //printf("-------------------------\nBoundary Element readed correctly \n-------------------------\n");    
    // LOAD GEOM IDS
    //--------------------------
    std::getline(inFile, line);
    std::getline(inFile, line); std::getline(inFile, line);
    std::shared_ptr<int[]> ids(new int[NBoundelem]);
    int incumb_max_id = 0;

    //printf("-------------------------\nStarting Boundary ID read \n-------------------------\n");    
    for(int i = 0; i < NBoundelem; i++)
    {
        std::getline(inFile, line);
        iss.clear();
        iss.str(line);
        int curr_id;
        iss >> curr_id; ids[i] = curr_id;
        if (curr_id > incumb_max_id) incumb_max_id = curr_id;
    }
    //printf("-------------------------\nBoundary ID readed correctly\n-------------------------\n");   
    int nBound = incumb_max_id+1;
    std::vector<VECTOR_INT> ElemOrderedById(nBound);
    VECTOR_INT idCount(nBound); idCount.zeros();

    // ALLOC FOR EACH BOUND
    for (int i = 0; i < nBound; i++) {
        ElemOrderedById[i].initialize(NBoundelem);  
    }
    for (int iel = 0; iel < NBoundelem; iel++)
    {
        int id = ids[iel];
        int locCount = idCount[id];
        ElemOrderedById[id][locCount] = iel; 
        idCount[id] ++;
    }
    // print Bound Elems by ID to file
    BoundElemFile << NBoundelem << " " << nBound << "\n";
    for (int id = 0; id < nBound; id++)
    {
        BoundElemFile << id << " " << idCount[id] << "\n";
        for (int iel = 0; iel < idCount[id]; iel++)
        {
            int elem = ElemOrderedById[id][iel];
            int* tempNodes = BoundElems[elem];
            for (int inode = 0; inode < boundNodesxElem; inode++)
            {
                int locNode =  tempNodes[inode];
                BoundElemFile << locNode << " ";
            }
            BoundElemFile << "\n";
        }
        BoundElemFile << "\n";

        ElemOrderedById[id].shrink(idCount[id]);
    }
    BoundElemFile.close();
    //-----------------------------------
    // PRINT BOUNDARY NODES BY ID TO FILE
    //-----------------------------------
    VECTOR_INT BoundNodesByID(NBoundelem*boundNodesxElem);
    int globCount = 0;
    VECTOR_INT locCount(nBound);
    for (int id = 0; id < nBound; id++)
    {
        VECTOR_INT isVisited(Nodes); // for bound nodes only
        int locCountNodes = 0;
        for (int iel = 0; iel < idCount[id]; iel++)
        {
            int elem = ElemOrderedById[id][iel];
            int* tempNodes = BoundElems[elem];
            for (int inode = 0; inode < boundNodesxElem; inode++)
            {
                int locNode =  tempNodes[inode];
                isVisited[locNode] = 0;
            }
        }
        //---
        for (int iel = 0; iel < idCount[id]; iel++)
        {
            int elem = ElemOrderedById[id][iel];
            int* tempNodes = BoundElems[elem];
            for (int inode = 0; inode < boundNodesxElem; inode++)
            {
                int locNode =  tempNodes[inode];
                if (isVisited[locNode] == 0)
                {
                    isVisited[locNode] = 1; // set to visited
                    BoundNodesByID[globCount] = locNode;
                    locCountNodes++;
                    globCount++;
                }
            }
        }
        locCount[id] = locCountNodes;
    }
    // print values to file
    std::vector<VECTOR_INT> boundNode(nBound);
    BoundNodeFile << globCount << " " << nBound << "\n";
    int count = 0;
    for (int id = 0; id < nBound; id++)
    {
        int locCountNodes = locCount[id];
        boundNode[id].initialize(locCountNodes);
        std::shared_ptr<int[]> locPointer = boundNode[id].P;

        BoundNodeFile << id << " " << locCountNodes << "\n";
        for (int inode = 0; inode < locCountNodes; inode++)
        {
            int tempNode = BoundNodesByID[count];
            BoundNodeFile << tempNode << "\n";
            locPointer[inode] = tempNode;
            count++;
        }
        BoundNodeFile << "\n";
    }

    //-------------
    // PRINT ON VTK
    //-------------
    VTK vtkWriter;
    vtkWriter.initFileBounds(name, nDim, nBound);
    int countBNode = 0;
    for (int ibound = 0; ibound < nBound; ibound++)
    {
        int nNode = locCount[ibound]; 
        MATRIX tempCoord(nNode, nDim);
        //---
        VECTOR_INT locBoundNodes(nNode);
        //---
        for (int j = 0; j < nNode; j++)
        {
            int locNode = BoundNodesByID[countBNode];
            prec* tempNode = tempCoord[j];

            prec* tempP = coord[locNode];
            VECTOR::copy(tempP, tempNode, nDim);
            locBoundNodes[j] = locNode;
            countBNode++;
        }
        //---
        int tempNelem = idCount[ibound];
        MATRIX_INT tempElems(tempNelem, nDim);
        std::shared_ptr<int[]> localStartElem = ElemOrderedById[ibound].P;
        for (int j = 0; j < tempNelem; j++)
        {
            int locElem = localStartElem[j];
            int* tempNodes = BoundElems[locElem];
            int* tempElem = tempElems[j];
            for (int k = 0; k < nDim; k++)
            {
                int tempNode = tempNodes[k];
                int index;
                bool found = locBoundNodes.hasIn(tempNode, index);
                if (!found) throw_line("ERROR: wrong writing of BoundNodes\n");
                tempElem[k] = index;
            }
        }
        vtkWriter.writeBound(ibound, tempCoord, tempElems);
    }
    vtkWriter.closeBoundFile();

    // skip lines
    //--------------------------
    for (int i = 0; i < 6; i++) std::getline(inFile, line);


    //--------------------------
    // LOAD ELEMS
    //-------------------------
    std::getline(inFile, line);
    iss.clear();
    iss.str(line);
    //---------------------------
    int NodesxElem;
    iss >> NodesxElem;
    int Nelem;
    std::getline(inFile, line);
    iss.clear();
    iss.str(line); iss >> Nelem;
    std::getline(inFile, line);
    ElemFile << Nelem << " \n";
    // countPerc = 1; maxCountPerc = Nelem/1000;

    MATRIX_INT elem(Nelem, NodesxElem);
    //printf("-------------------------\nStarting Elements load \n-------------------------\n");   
    for(int i = 0; i < Nelem; i++)
    {
        std::getline(inFile, line);
        iss.clear();
        iss.str(line);
        int* tempElem = elem[i];
        for (int j = 0; j < NodesxElem; j++)
        {
            int value;
            iss >> value;
            ElemFile << value << " ";

            tempElem[j] = value;
        }
        ElemFile << "\n ";

        //print percentage of loading
        // if (countPerc == maxCountPerc)
        // {
        //     prec tempPerc = (prec) 100*(i+1)/Nelem;
        //     //printf("Loading Elements: %4.2" format "%% \n", tempPerc);
        //     countPerc = 0;
        // }
        // countPerc++;
    }
    //printf("-------------------------\nElements loaded correctly \n-------------------------\n");   
    ElemFile.close();

    //------READ GEOMETRIC ENETITIES INDICES------
    // This is done in order to be able to define the initial value of gamma in Top_Opt problems in different domain regions
    std::getline(inFile, line);
    int n_geo_entities_ids;
    STREAM::getValue(inFile, line, iss, n_geo_entities_ids);
    if (n_geo_entities_ids != Nelem) throw_line("ERROR: different number of elements and geometric entity indices.\n");
    std::getline(inFile, line);
    VECTOR_INT geo_entities_ids(Nelem);
    STREAM::getColVector(inFile, line, iss, geo_entities_ids, Nelem);
    for (int iel = 0; iel < Nelem; iel++)
    {
        geo_entities_ids[iel] -= 1; // use C++ notation starting from 0
    }

    // print on file
    ElemGeoIdsFile << Nelem << "\n";
    for (int i = 0; i < Nelem; i++)
    {
        ElemGeoIdsFile << geo_entities_ids[i] << "\n";
    }
    ElemGeoIdsFile.close();
    //printf("-------------------------\nElements indices loaded correctly \n-------------------------\n");  

    inFile.close();
    free(folderName);




    //----------------------------------------
    // CREATE P1-ISO P2/P1 MESH
    //----------------------------------------
    //printf("11\n");
    int Nelem_v = Nelem*(4*nDim-4);
    int nEvals = (4*nDim-2)*Nelem;
    std::shared_ptr<int[]> i_sparse(new int[nEvals]);
    std::shared_ptr<int[]> j_sparse(new int[nEvals]);
    std::shared_ptr<prec[]> coef_sparse(new prec[nEvals]);

    int countForAss = 0;
    for (int iel = 0; iel < Nelem; iel++)
    {
        int* tempEl = elem[iel];
        for (int iloc = 0; iloc < nDim; iloc++)
        {
            int iglob = tempEl[iloc];
            for (int jloc = iloc; jloc < nDim+1; jloc++)
            {
                int jglob = tempEl[jloc];
                if (iglob > jglob)
                {
                    // if (iglob == 2)
                    // {
                    //     printf("iglob: %d, jglob: %d\n", iglob, jglob);
                    // }
                    i_sparse[countForAss] = iglob;
                    j_sparse[countForAss] = jglob;
                }
                else
                {
                    // if (jglob == 2)
                    // {
                    //     printf("iglob: %d, jglob: %d\n", jglob, iglob);
                    // }
                    i_sparse[countForAss] = jglob;
                    j_sparse[countForAss] = iglob;
                }
                coef_sparse[countForAss] = 0;
                countForAss++;
            }
        }
    }
    CSRMAT isCreated(Nodes, Nodes, countForAss, i_sparse, j_sparse, coef_sparse);
    //isCreated.writePatt("this.txt");
    //printf("11\n");

    MATRIX coord_v(Nodes+Nelem*3*(nDim-1), nDim);
    MATRIX_INT elem_v(Nelem_v, nDim+1);
    VECTOR_INT geo_entities_ids_v(Nelem_v);
    

    std::vector<MATRIX_INT> boundElemById_v(nBound);
    // alloc space for bound velocity elements
    int factor = 2*(nDim-1);
    for (int ibound = 0; ibound < nBound; ibound++)
    {
        boundElemById_v[ibound].initialize(factor*idCount[ibound], nDim);
    }
    VECTOR_INT boundElCounter = VECTOR_INT::zeros(nBound);
    int currNod = Nodes;

    for (int i = 0; i < Nodes; i++)
    {
        prec* tempCoord = coord_v[i];
        prec * tempP = coord[i];
        VECTOR::copy(tempP, tempCoord, nDim);
    } 
    //---------------------------------ù
    // ALLOC BOUND NODES
    //---------------------------------
    //printf("0\n");
    std::vector<VECTOR_INT> boundNode_v(nBound);
    std::vector<int> boundNodesCounter(nBound);
    for (int ibound = 0; ibound < nBound; ibound++)
    {
        int oldLength = boundNode[ibound].length;
        int newLength = oldLength+idCount[ibound]*3;
        
        boundNode_v[ibound].initialize(newLength);
        std::shared_ptr<int[]> tempOldPointer = boundNode[ibound].P;
        std::shared_ptr<int[]> tempNewPointer = boundNode_v[ibound].P;
        for (int inode = 0; inode < oldLength; inode++)
        {
            tempNewPointer[inode] = tempOldPointer[inode];
        }

        boundNodesCounter[ibound] = oldLength;
    }
    /// start from bound elements
    /// BOUND ELEMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //---------------------------------------------------------------------------------------
    int locLimit = 2*nDim-3;


    CSRMAT tempCreated;
    isCreated.copyPattTo(tempCreated);
    

    // LOWER_TRIANGULAR tempCreated(Nodes);

    for (int id = 0; id < nBound; id++)
    {
        int tempCounter = 0;
        // printf("boundID: %d\n", id);
        for (int i = 0; i < tempCreated.nTerm; i++) tempCreated.coef[i] = 0;
        
        for (int iel = 0; iel < idCount[id]; iel++)
        {
            int elem = ElemOrderedById[id][iel];
            int* vertices = BoundElems[elem];
            std::vector<VECTOR> tempCoords(boundNodesxElem);
            std::vector<int> tempId(locLimit);
            for (int iloc = 0; iloc < boundNodesxElem; iloc++)
            {
                tempCoords[iloc].initialize(nDim);
                int iglob = vertices[iloc];
                prec* tempP = coord[iglob];
                tempCoords[iloc] = tempP;
            }
            //---
            /// tempId(1) = middle of arc 1-2
            /// tempId(2) = middle of arc 2-3
            /// tempId(3) = middle of arc 3-1
            //---
            for (int iloc = 0; iloc < locLimit; iloc++)
            {
                int jloc = (iloc+1)%boundNodesxElem;
                int iglob = vertices[iloc]; 
                int jglob = vertices[jloc]; 
                prec* created;
                int posGlob;
                // printf("IGLOB: %d, JGLOB: %d\n", iglob, jglob);
                if (iglob > jglob)
                {
                    posGlob = isCreated.getPos(iglob, jglob);
                    created = &(isCreated.coef[posGlob]);
                }
                else 
                {
                    posGlob = isCreated.getPos(jglob, iglob);
                    created = &(isCreated.coef[posGlob]);
                }
                // CHECK IF CREATED
                // if not created
                if (*created == 0)
                {
                    VECTOR newPoint = (tempCoords[iloc]+tempCoords[jloc])/2;
                    *created = currNod;
                    prec* temp = coord_v[currNod];
                    newPoint.copyTo(temp);
                    tempId[iloc] = currNod;

                    boundNode_v[id][boundNodesCounter[id]] = currNod;
                    boundNodesCounter[id]++;
                    // create boundary elements with appropriate bound id

                    if (iglob > jglob) tempCreated.coef[posGlob] = 1;
                    else tempCreated.coef[posGlob] = 1;
                    currNod++;
                } else
                //--
                // if created
                {
                    tempId[iloc] = *created;
                    if (iglob > jglob && tempCreated.coef[posGlob] != 1) 
                    {
                        tempCreated.coef[posGlob] = 1;
                        boundNode_v[id][boundNodesCounter[id]] = *created;
                        boundNodesCounter[id]++;
                    }
                    else if (jglob > iglob && tempCreated.coef[posGlob] != 1 )
                    {
                        tempCreated.coef[posGlob] = 1;
                        boundNode_v[id][boundNodesCounter[id]] = *created;
                        boundNodesCounter[id]++;
                    } 
                }
            }
            // CREATE BOUNDARY ELEMENTS
            if (nDim == 2)
            {
                std::shared_ptr<int*[]> tempBoundElem = boundElemById_v[id].PP;
                int* tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set2Values(tempElem, vertices[0], tempId[0]);
                tempCounter++;
                tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set2Values(tempElem, tempId[0], vertices[1]);
                tempCounter++;
            } else
            {
                std::shared_ptr<int*[]> tempBoundElem = boundElemById_v[id].PP;
                int* tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set3Values(tempElem, vertices[0], tempId[0], tempId[2]);
                tempCounter++;
                tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set3Values(tempElem, vertices[1], tempId[1], tempId[0]);
                tempCounter++;
                tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set3Values(tempElem, vertices[2], tempId[2], tempId[1]);

                tempCounter++;
                tempElem = tempBoundElem[tempCounter];
                VECTOR_INT::set3Values(tempElem,  tempId[0], tempId[1], tempId[2]);
                tempCounter++;
            }
        }
        boundElCounter[id] = tempCounter;
    }
    // printf("3\n");
    tempCreated.coef.reset(); tempCreated.ja.reset();

    // CREATE REMAINING MESH OF P1-ISO P2/P1
    locLimit = 3*nDim - 3;
    int elCounter = 0;
    for(int iel = 0; iel < Nelem; iel++)
    {
        //----
        int* vertices = elem[iel];
        std::vector<VECTOR> tempCoords(NodesxElem);
        std::vector<int> tempId(locLimit);
        for (int iloc = 0; iloc < NodesxElem; iloc++)
        {
            tempCoords[iloc].initialize(nDim);
            int iglob = vertices[iloc];
            prec* tempP = coord[iglob];
            tempCoords[iloc] = tempP;
        }
        //---
        /// tempId(1) = middle of arc 1-3
        /// tempId(2) = middle of arc 1-4
        /// tempId(3) = middle of arc 2-4
        /// tempId(4) = middle of arc 2-1
        /// tempId(5) = middle of arc 3-1
        /// tempId(6) = middle of arc 3-4
        //---
        
        // CASE 2 DIM
        if (nDim == 2)
        {
            for (int iloc = 0; iloc < 3; iloc++)
            {
                int jloc = (iloc+1)%NodesxElem;
                int iglob = vertices[iloc]; 
                int jglob = vertices[jloc]; 
                prec* created;
                int posGlob;
                if (iglob > jglob)
                {
                    posGlob = isCreated.getPos(iglob, jglob);
                    created = &(isCreated.coef[posGlob]);
                }
                else 
                {
                    posGlob = isCreated.getPos(jglob, iglob);
                    created = &(isCreated.coef[posGlob]);
                }
                // CHECK IF CREATED
                // if not created
                if (*created == 0)
                {
                    VECTOR newPoint = (tempCoords[iloc]+tempCoords[jloc])/2;
                    *created = currNod;
                    prec* temp = coord_v[currNod];
                    newPoint.copyTo(temp);
                    tempId[iloc] = currNod;
                    // create boundary elements with appropriate bound id
                    currNod++;
                } else
                {
                    tempId[iloc] = *created;
                }
            }
        } else // CASE 3D
        {
            int locCount = 0;
            for (int iloc = 0; iloc < NodesxElem; iloc++)
            {
                int iglob = vertices[iloc];
                for (int jloc = iloc+1; jloc < NodesxElem; jloc++)
                {
                    int jglob = vertices[jloc]; 
                    prec* created;
                    int posGlob;
                    if (iglob > jglob)
                    {
                        posGlob = isCreated.getPos(iglob, jglob);
                        created = &(isCreated.coef[posGlob]);
                    }
                    else 
                    {
                        posGlob = isCreated.getPos(jglob, iglob);
                        created = &(isCreated.coef[posGlob]);
                    }
                    // CHECK IF CREATED
                    // if not created
                    if (*created == 0)
                    {
                        VECTOR newPoint = (tempCoords[iloc]+tempCoords[jloc])/2;
                        *created = currNod;
                        prec* temp = coord_v[currNod];
                        newPoint.copyTo(temp);
                        tempId[locCount] = currNod;
                        // create boundary elements with appropriate bound id
                        currNod++;
                    } else // if created
                    {
                        tempId[locCount] = *created;
                    }
                    locCount++;
                }
            }
        }

        // CREATE ELEMENTS
        int temp_elem_geo_id = geo_entities_ids[iel];
        if (nDim == 2)
        {
            int* tempElem = elem_v[elCounter];
            VECTOR_INT::set3Values(tempElem, vertices[0], tempId[0], tempId[2]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set3Values(tempElem, vertices[1], tempId[1], tempId[0]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set3Values(tempElem, vertices[2], tempId[2], tempId[1]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set3Values(tempElem,  tempId[0], tempId[1], tempId[2]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
        } else
        {
            int* tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, vertices[0], tempId[0], tempId[1], tempId[2]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, tempId[0], vertices[1], tempId[3], tempId[4]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, tempId[1], tempId[3], vertices[2], tempId[5]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, tempId[2], tempId[4], tempId[5], vertices[3]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, tempId[2], tempId[1], tempId[4], tempId[0]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter]; 
            VECTOR_INT::set4Values(tempElem, tempId[1], tempId[3], tempId[4], tempId[0]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter];
            VECTOR_INT::set4Values(tempElem, tempId[5], tempId[4], tempId[1], tempId[3]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
            tempElem = elem_v[elCounter]; 
            VECTOR_INT::set4Values(tempElem, tempId[5], tempId[4], tempId[2], tempId[1]);
            geo_entities_ids_v[elCounter] = temp_elem_geo_id;
            elCounter++;
        }
        //-------------
    }
    std::string pathCoord_v = folderPath + "/Nodes_v.txt";
    coord_v.print3(&pathCoord_v[0], currNod);
    std::string pathElem_v = folderPath + "/Elems_v.txt";
    elem_v.print3(&pathElem_v[0]);
    std::string pathElemGeoIds_v = folderPath + "/ElemGeoIds_v.txt";
    geo_entities_ids_v.print(&pathElemGeoIds_v[0], true);

    std::string boundFileName = folderPath + "/BoundNodes_v.txt";
    std::ofstream boundFile_v; boundFile_v.open(boundFileName);

    int sum = 0;
    for (int id = 0; id < nBound; id++)
    {
        sum += boundNodesCounter[id];

    }
    boundFile_v << sum << " " << nBound << "\n";
    for (int id = 0; id < nBound; id++)
    {
        int tempLength = boundNodesCounter[id];
        std::shared_ptr<int[]> tempPointer = boundNode_v[id].P;
        boundFile_v << id << " " << tempLength << "\n";
        for (int i = 0; i < tempLength; i++)
        {
            boundFile_v << tempPointer[i] << "\n";
        }
        boundFile_v << "\n";
    }
    boundFile_v.close();
    //-------------------------------
    // PRINT BOUND ELEMS
    //-------------------------------
    ////printf("PRINT BOUND ELEMS \n");
    std::string boundElemFileName = folderPath + "/BoundElems_v.txt";
    std::ofstream boundElemFile_v; boundElemFile_v.open(boundElemFileName);
    sum = 0;
    for (int id = 0; id < nBound; id++)
    {
        sum += boundElCounter[id];

    }
    boundElemFile_v << sum << " " << nBound << "\n";
    for (int id = 0; id < nBound; id++)
    {
        int tempLength = boundElCounter[id];
        boundElemFile_v << id << " " << tempLength << "\n";
        //----------------------------------------------
        int** tempMat = &(boundElemById_v[id].PP[0]);
        for (int iel = 0; iel < tempLength; iel++)
        {
            int* tempCoord = tempMat[iel];
            for (int j = 0; j < nDim; j++)
            {
                boundElemFile_v << tempCoord[j] << " ";
            }
            boundElemFile_v << "\n";
        }
        boundElemFile_v << "\n";
    }
    boundElemFile_v.close();

    int nDof = currNod*nDim + Nodes;
    return nDof;
}

