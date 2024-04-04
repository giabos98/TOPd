#pragma once

#include "codeHeader.h"
namespace fs = std::filesystem;

class VTK
{
public:
    std::string path;
    std::string currFilePath;
    std::string problemName;
    bool binWrite = false;
    int flagPrint;
    int dim;
    int nNode;
    int nElem;
    int nPrint;
    int counter = 0;
    int counterPDE = 0;
    int max_opt_it = 0;

    prec currPrintTime;
    prec deltaT;
    std::shared_ptr<prec[]> exactTimes;

    std::ofstream boundFile;
    std::ofstream tFile;
    std::string path_Tseries_file;
    std::string tString = "";
    std::string close_tString = "\n\t]\n}";
    std::string frmt;
    //------------------------
    // CONSTRUCTOR
    //------------------------
    VTK(){}
    //-------------
    void initialize(int dim_in, int nNodes, int nElems, std::string name, bool bin = false)
    {
        problemName = name;
        dim = dim_in;
        nNode = nNodes;
        nElem = nElems;

        path = "../testBinary/" + name;

        binWrite = bin;

        if (binWrite) frmt = "binary";
        else frmt = "ascii";
    }
    //-------------
    void initializeForTopOpt(std::string name, int dim_in, int nNodes, int nElems, int deltaPrint, int maxIt, bool bin)
    {
        nNode = nNodes; nElem = nElems; dim = dim_in;
        path = "results/" + name;
        fs::create_directories(&path[0]);

        problemName = name;
        // binread
        binWrite = bin;

        //------------
        if (binWrite) frmt = "binary";
        else frmt = "ascii";

        deltaT = deltaPrint;
        max_opt_it = maxIt;

        openTSeriesFile();
        currPrintTime = 1;
    }
    //----------------------------------------------------------------------
    void initializeForNS(std::string name, int dim_in, int nNodes, int nElems, std::string folderName = "NS_sol")
    {
        nNode = nNodes; nElem = nElems; dim = dim_in;
        name = name + "/" + folderName;
        path = "results/" + name;
        fs::create_directories(&path[0]);

        problemName = name;
        //----- read print infos from "infoPrintResults.txt" ----------------
        std::string infoPath  = "INPUT_FILES/infoPrintNS.txt";
        //-------------------------
        std::ifstream  file;  file.open(infoPath, std::ios::in);  if (!file.is_open()) throw_line("ERROR, can't open infoPrintResults file");
        std::string line;
        std::istringstream iss;
        // binread
        std::getline(file, line);
        std::getline(file, line);
        iss.str(line);
        iss >> binWrite;
        iss.clear();
        // nPrint
        STREAM::getLines(file, line, 2);
        STREAM::getValue(file, line, iss, nPrint);
        // flagPrint
        STREAM::getLines(file, line, 4);
        STREAM::getValue(file, line, iss, flagPrint);
        // values

        if (flagPrint == 0) // read t start and  delta T 
        {
            STREAM::getLines(file, line, 2);
            STREAM::getValue(file, line, iss, currPrintTime);
            STREAM::getLines(file, line, 1);
            STREAM::getValue(file, line, iss, deltaT);
        } 
        else
        {
            exactTimes = std::shared_ptr<prec[]>(new prec[nPrint+1]);
            exactTimes[0] = 0;
            prec* tempPoint = &exactTimes[1];
            STREAM::getLines(file, line, 7);
            STREAM::getColVector(file, line, iss, tempPoint, nPrint);
            currPrintTime = exactTimes[0];
        }
        file.close();
        //------------
        if (binWrite) frmt = "binary";
        else frmt = "ascii";
    }
    //------------------------------------------------------------------------
    void initializeForADJ(std::string name, int dim_in, int nNodes, int nElems, std::string folderName = "ADJ_sol")
    {
        initializeForNS(name, dim_in, nNodes, nElems, folderName);
    }
    //-----------------
    // RESET FOR NS
    //-----------------
    void resetForNS(std::string newFolderName)
    {
        path = "results/" + problemName + "/" + newFolderName;
        fs::create_directories(&path[0]);

        currPrintTime = 0;
        openTSeriesFile();
    }
    //---------------------------------
    // PRINT VELOCITIES?
    void writeMesh(int nNode, int nElem, MATRIX &coord, MATRIX_INT &elems)
    {
        std::string currFileName = "mesh.vtu";
        currFilePath = path + "/" + currFileName;
        std::ofstream file;
        if (binWrite) file.open(currFilePath, std::ios::out | std::ios::binary);
        else file.open(currFilePath);

        file << "<?xml version=\"1.0\"?>\n";
        //--- UNSTRUCTURED GRID ---------------
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        file << "  <UnstructuredGrid>\n";
        file << "\t\t<Piece NumberOfPoints=\"" << nNode << "\" NumberOfCells=\"" << nElem << "\">\n";

        writeMesh(file, coord, elems);
        file << "\t\t</Piece>\n";

        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>";
    }
    //---------------------------------
    template <class... Types>
    void write(MATRIX &coord, MATRIX_INT &elems, prec time, Types... args)
    {
        int curr_it = counter+1;
        int delta_len = std::to_string(max_opt_it).length() - std::to_string(curr_it).length();
        std::string count_str = "";
        for (int ichar = 0; ichar < delta_len; ichar++)
        {
            count_str += "0";
        }
        count_str += std::to_string(curr_it);
        std::string currFileName = "it" + count_str + ".vtu";
        currFilePath = path + "/" + currFileName;
        std::ofstream file;
        if (binWrite) file.open(currFilePath, std::ios::out | std::ios::binary);
        else file.open(currFilePath);

        file << "<?xml version=\"1.0\"?>\n";
        //--- UNSTRUCTURED GRID ---------------
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
        file << "  <UnstructuredGrid>\n";
        writeTime(file, time);
        file << "\t\t<Piece NumberOfPoints=\"" << nNode << "\" NumberOfCells=\"" << nElem << "\">\n";
        //-----------------------------------------
        int nArgs = sizeof...(args);
        if (nArgs%3 != 0) throw_line("Missing parameter in VTK writer\n");
        writeSol(file, args...);

        writeMesh(file, coord, elems);
        file << "\t\t</Piece>\n";

        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>";

        // if (counter != 0) 
        // {
        //     tFile << ",\n";
        //     tString += ",\n";
        // }
        TSeriesWriteTime(currFileName, time);

        // counter++;
    }
    //---
    template <class... Types>
    void writeTimeDep(MATRIX &coord, MATRIX_INT &elems, prec time, prec currDeltaT, Types &...args)
    {
        if (counter > nPrint) return;
        if (abs(time - currPrintTime) > abs(currPrintTime - time - currDeltaT)) return;
        write(coord, elems, time, args...);
        currPrintTime += currDeltaT;
    }
    //------------------------------------------------
    void closeTFile()
    {
        // tString += close_tString;
        // tFile << "\n\t]\n";
        // tFile << "}";
        // tFile.close();
    }
    //----------------------------------------------------
    // BOUND ----------------------- BOUND----------------
    //----------------------------------------------------
    void initFileBounds(std::string name, int dim_in, int nBound)
    {
        std::string filePath = "PROBLEM_DATA/" + name;
        frmt = "ascii";
        std::string pathDir = filePath + "/plotBound";
        filePath = pathDir + "/bound.pvtu";

        fs::create_directories(&pathDir[0]);
        binWrite = false;
        if (binWrite) boundFile.open(filePath, std::ios::binary);
        else  boundFile.open(filePath);

        dim = dim_in;
        path = pathDir;

        boundFile << "<?xml version=\"1.0\"?>\n";
        //--- UNSTRUCTURED GRID ---------------
        boundFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        boundFile << "<PUnstructuredGrid GhostLevel=\"" << nBound << "\">\n";
    }
    //-----------
    void writeBound(int idBound, MATRIX &coord, MATRIX_INT &elem)
    {
        int nPoints = coord.nRow; int nElem = elem.nRow;
        std::string tempName = "bound" + std::to_string(idBound) +".vtu";
        boundFile << "<PPoints>\n";
        boundFile << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
        boundFile << "</PPoints>\n";
        boundFile << "<Piece Source=\"" << tempName << "\"/>\n";

        std::ofstream tempFile;
        std::string tempFilePath = path + "/" + tempName;
        if (binWrite) tempFile.open(tempFilePath, std::ios::binary);
        else  tempFile.open(tempFilePath);

        tempFile << "<?xml version=\"1.0\"?>\n";
        //--- UNSTRUCTURED GRID ---------------
        tempFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        tempFile << " <UnstructuredGrid>\n";

        tempFile << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nElem << "\">\n";
        writeMesh(tempFile, coord, elem);
        tempFile << "\t\t</Piece>\n";

        tempFile << "  </UnstructuredGrid>\n";
        tempFile << "</VTKFile>";
        tempFile.close();
    }
    //--------
    void closeBoundFile()
    {
        boundFile << "</PUnstructuredGrid>\n";
        boundFile << "</VTKFile>";
        boundFile.close();
    }
private:
    //-----------
    void writeTime(std::ofstream &file, prec time)
    {
        file << "\t<FieldData>\n";
        file << "\t\t<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"" << frmt << "\" RangeMin=\"" << time << "\" RangeMax=\" " << time  <<"\">\n";
        file << "\t\t\t" << time << "\n";
        file << "\t\t</DataArray>\n";
        file << "\t</FieldData>\n";
    }
    //------------------------
    // write SOL and VELOCITY
    //------------------------
    void writeSolVel(std::string fileName, std::ofstream &file, VECTOR &sol, MATRIX &nodalVel)
    {
        file << "\t\t\t<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n";
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"" << frmt << "\">\n";
        file.close();

        FILE* outFile = fopen(&fileName[0], "a");
        for (int i = 0; i < nNode; i++)
        {
            fprintf(outFile, "\t\t\t\t\t");
            prec* tempVel = nodalVel[i];
            if (dim == 2)
            {
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    prec tempVal = tempVel[icomp];
                    if (abs(tempVal) < 1e-17) fprintf(outFile, "0 "); //file << 0 << " ";
                    else fprintf(outFile, "%" format_e " ", tempVal); //file << tempVal << " ";
                }
                fprintf(outFile, "0\n");
            }
            else
            {
                for (int icomp = 0; icomp < dim-1; icomp++)
                {
                    prec tempVal = tempVel[icomp];
                    if (abs(tempVal) < 1e-17) fprintf(outFile, "0 "); //file << 0 << " ";
                    else fprintf(outFile, "%" format_e " ", tempVal); //file << tempVal << " ";
                }
                prec tempVal = tempVel[dim-1];
                if (abs(tempVal) < 1e-17) fprintf(outFile, "0\n"); //file << 0 << " ";
                else fprintf(outFile, "%" format_e "\n", tempVal); //file << tempVal << " ";
            }
        }
        fclose(outFile);
        file.open(fileName, std::ios_base::app);
        file << "\t\t\t\t</DataArray>\n";
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" format=\"" << frmt << "\">\n";

        for (int i = 0; i < nNode; i++)
        {
            file << "\t\t\t\t\t";
            prec tempSol = sol[i];
            if (abs(tempSol) < 1e-17) file << 0 << "\n";
            else file << tempSol << "\n";
        }
        file << "\t\t\t\t</DataArray>\n";
        file << "\t\t\t</PointData>\n";
        file << "\t\t\t<CellData>\n\t\t\t</CellData>\n";
    }
    //-- binary
    void writeSolVelBin(std::ofstream &file, VECTOR &sol, MATRIX &nodalVel)
    {
        file << "\t\t\t<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n";
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"binary\">\n";

        // std::ofstream newFile; newFile.open("prova.txt", std::ios::out | std::ios::binary);
        for (int i = 0; i < nNode; i++)
        {
            prec* tempVel = nodalVel[i];
            if (dim == 2)
            {
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    float tempVal = (float) tempVel[icomp];
                    if (abs(tempVal) < 1e-17) tempVal = 0;
                    file.write((char*)&tempVal, sizeof(float));
                }
                float pp = 0;
                file.write((char*)&pp, sizeof(float));
            }
            else
            {
                for (int icomp = 0; icomp < dim-1; icomp++)
                {
                    prec tempVal = tempVel[icomp];
                    if (abs(tempVal) < 1e-17) file << 0 << " ";
                    else file << tempVal << " ";
                }
                prec tempVal = tempVel[dim-1];
                if (abs(tempVal) < 1e-17) file << 0 << "\n";
                else file << tempVal << "\n";
            }
        }
        // newFile.close();
        // std::ifstream oldFile; oldFile.open("prova.txt", std::ios::in | std::ios::binary);
        // for (int i = 0; i < nNode; i++)
        // {
        //     if (dim == 2)
        //     {
        //         for (int icomp = 0; icomp < 3; icomp++)
        //         {
        //             float tempVal;
        //             oldFile.read(reinterpret_cast<char*>(&tempVal), sizeof(float));
        //             std::cout << tempVal << " ";
        //         }
        //         std::cout << "\n";
        //     }
        // }
        // oldFile.close();



            //------------

        file << "\t\t\t\t</DataArray>\n";
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" format=\"" << frmt << "\" RangeMin=\"0\" RangeMax=\""<< 20 << "\">\n";

        for (int i = 0; i < nNode; i++)
        {
            file << "\t\t\t\t\t";
            prec tempSol = sol[i];
            if (abs(tempSol) < 1e-17) file << 0 << "\n";
            else file << tempSol << "\n";
        }
        file << "\t\t\t\t</DataArray>\n";
        file << "\t\t\t</PointData>\n";
        file << "\t\t\t<CellData>\n\t\t\t</CellData>\n";
    }
    //-----------------------------------------------------
    template <class... Types>
    void writeSol(std::ofstream &file, Types&... args)
    {
        std::string ScalarText = "\"";
        std::string VectorText = "\"";
        bool foundVector = false;
        bool foundScalar = false;

        writeSolHeader(ScalarText, VectorText, foundScalar, foundVector, args...);

        ScalarText += "\"";
        VectorText += "\"";
        file << "\t\t\t<PointData Scalars=" + ScalarText + " Vectors=" + VectorText + ">\n";

        writeSol2(file, args...);

        file << "\t\t\t</PointData>\n";
        file << "\t\t\t<CellData>\n\t\t\t</CellData>\n";
    }
    //----------
    template <typename T, class... Types>
    void writeSolHeader(std::string &ScalarText, std::string & VectorText, bool &foundScalar, bool &foundVector, T &sol, int dim, std::string name, Types&... args)
    {

        switch (dim)
        {
            case 1:
            {
                if (foundScalar) ScalarText += ", ";
                ScalarText += name;
                foundScalar = true;
                break;
            }
            case 2:
            {
                if (foundVector) VectorText += ", ";
                VectorText += name;
                foundVector = true;
                break;
            }
            case 3:
            {
                if (foundVector) VectorText += ", ";
                VectorText += name;
                foundVector = true;
                break;
            }
            default:
            {
                std::string errMsg = "ERROR: wrong parameter " + name + " dimensions in VTK writer input\n";
                throw_line(&errMsg[0]);
            }
        }
        writeSolHeader(ScalarText, VectorText, foundScalar, foundVector, args...);
    }
    void writeSolHeader(std::string &ScalarText, std::string & VectorText, bool &foundScalar, bool &foundVector) {}
    //------------------------
    template <class... Types>
    void writeSol2(std::ofstream &file, VECTOR &sol, int dim, std::string name, Types&... args)
    {
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"" << name << "\" format=\"" << frmt << "\">\n";

        for (int i = 0; i < nNode; i++)
        {
            file << "\t\t\t\t\t";
            prec tempSol = sol[i];
            if (abs(tempSol) < 1e-17) file << 0 << "\n";
            else file << tempSol << "\n";
        }
        file << "\t\t\t\t</DataArray>\n";
        writeSol2(file, args...);
    }
    //---
    template <class... Types>
    void writeSol2(std::ofstream &file, MATRIX &sol, int dim, std::string name, Types&... args)
    {
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"" << frmt << "\">\n";
        file.close();
        FILE* outFile = fopen(&currFilePath[0], "a");
        for (int i = 0; i < nNode; i++)
        {
            fprintf(outFile, "\t\t\t\t\t");
            prec* tempVel = sol[i];
            if (dim == 2)
            {
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    prec tempVal = tempVel[icomp];
                    if (abs(tempVal) < 1e-17) fprintf(outFile, "0 "); //file << 0 << " ";
                    else fprintf(outFile, "%" format_e " ", tempVal); //file << tempVal << " ";
                }
                fprintf(outFile, "0\n");
            }
            else
            {
                for (int icomp = 0; icomp < dim-1; icomp++)
                {
                    prec tempVal = tempVel[icomp];
                    if (abs(tempVal) < 1e-17) fprintf(outFile, "0 "); //file << 0 << " ";
                    else fprintf(outFile, "%" format_e " ", tempVal); //file << tempVal << " ";
                }
                prec tempVal = tempVel[dim-1];
                if (abs(tempVal) < 1e-17) fprintf(outFile, "0\n"); //file << 0 << " ";
                else fprintf(outFile, "%" format_e "\n", tempVal); //file << tempVal << " ";
            }
        }
        fclose(outFile);
        file.open(currFilePath, std::ios_base::app);
        file << "\t\t\t\t</DataArray>\n";
        writeSol2(file, args...);
    }
    //---
    void writeSol2(std::ofstream &file) {}
    //----------
    void writeMesh(std::ofstream &file, MATRIX &coord, MATRIX_INT &elem)
    {
        file << "\t\t\t<Points>\n";
        file << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << frmt << "\">\n";

        nNode = coord.nRow;
        if (dim == 2)
        {
            if (binWrite)
            {
                for (int i = 0; i < nNode; i++)
                {
                    prec* locCoord = coord[i];
                    for (int j = 0; j < dim; j++) 
                    {
                        float tempVar = locCoord[j];
                        file.write(reinterpret_cast<const char*>(&tempVar),sizeof(float));
                    }
                    float tempVar2 = 0;
                    file.write(reinterpret_cast<const char*>(&tempVar2),sizeof(float));
                }
                file << "\n";
            }
            else
            {
                for (int i = 0; i < nNode; i++)
                {
                    file << "\t\t\t\t\t";
                    prec* locCoord = coord[i];
                    for (int j = 0; j < dim; j++) file << locCoord[j] << " ";
                    file << 0; 
                    file << "\n";
                }
            }
        }
        else
        {
            if (binWrite)
            {
                for (int i = 0; i < nNode; i++)
                {
                    prec* locCoord = coord[i];
                    float tempVar;
                    for (int j = 0; j < dim; j++) 
                    {
                        tempVar = locCoord[j];
                        file.write(reinterpret_cast<const char*>(&tempVar),sizeof(float));
                    }
                    tempVar = locCoord[2];
                    file.write(reinterpret_cast<const char*>(&tempVar),sizeof(float));
                }
                file << "\n";
            }
            else
            {
                for (int i = 0; i < nNode; i++)
                {
                    file << "\t\t\t\t\t";
                    prec* locCoord = coord[i];
                    for (int j = 0; j < dim-1; j++) file << locCoord[j] << " ";
                    file << locCoord[2]; 
                    file << "\n";
                }
            }
        }
        file << "\t\t\t\t</DataArray>\n"; 
        file << "\t\t\t</Points>\n";
        //----------
        file << "\t\t\t<Cells>\n";
        
        nElem = elem.nRow; int nodesXElem = elem.nCol;
        file << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"" << frmt << "\" RangeMin=\"0\" RangeMax=\""<< nNode-1 << "\">\n";
        if (binWrite)
        {
            for (int i = 0; i < nElem; i++)
            {
            int* locNodes = elem[i];
            int tempVar;
            for (int j = 0; j < nodesXElem-1; j++)
            {
                tempVar = locNodes[j];
                file.write(reinterpret_cast<const char*>(&tempVar),sizeof(int));
            } 
            tempVar = locNodes[nodesXElem-1];
            file.write(reinterpret_cast<const char*>(&tempVar),sizeof(int));
            file << "\n";
            }
        }
        else
        {
            for (int i = 0; i < nElem; i++)
            {
            file << "\t\t\t\t\t";
            int* locNodes = elem[i];
            for (int j = 0; j < nodesXElem-1; j++) file << locNodes[j] << " ";
            file << locNodes[nodesXElem-1];
            file << "\n";
            }
        }
        
        file << "\t\t\t\t</DataArray>\n"; 
        //--
        int maxId = nElem*(nodesXElem);
        file << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"" << frmt << "\" RangeMin=\"" << nodesXElem << "\" RangeMax=\""<< maxId << "\">\n";
        int count = 0;
        if (binWrite)
        {
            for (int i = 0; i < nElem; i++)
            {
            count += nodesXElem;
            file.write(reinterpret_cast<const char*>(&count),sizeof(int));
            }
            file << "\n";
        }
        else
        {
            for (int i = 0; i < nElem; i++)
            {
            count += nodesXElem;
            file << "\t\t\t\t\t";
            file << count << "\n";
            }
        }
        
        file << "\t\t\t\t</DataArray>\n"; 
        //---
        int type = 3;
        switch (nodesXElem)
        {
        case 2: 
        {
            type = 3;
            break;
        }
        case 3:
        {
            type = 5;
            break;
        }
        case 4:
        {
            type = 10;
            break;
        }
        }
        file << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"" << frmt << "\" RangeMin=\"" << type << "\" RangeMax=\""<< type << "\">\n";
        if (binWrite)
        {
            for (int i = 0; i < nElem; i++)
            {
            file.write(reinterpret_cast<const char*>(&type),sizeof(int));
            }
            file << "\n";
        }
        else
        {
            for (int i = 0; i < nElem; i++)
            {
            file << type << "\n";
            }
        }
        
        file << "\t\t\t\t</DataArray>\n"; 
        //------
        file << "\t\t\t</Cells>\n";
        //----------
    }
    //----------------------------------------------------
    // TIME SERIES FILE (json file, extension .vtu.series)
    //----------------------------------------------------
    void openTSeriesFile()
    {
        path_Tseries_file = path + "/res.vtu.series";
        // if (binWrite) tFile.open(path_Tseries_file, std::ios::out | std::ios::binary);
        // else  tFile.open(path_Tseries_file);
        // tFile.close();
        tString = "";
        tString += "{\n\t\"file-series-version\" : \"1.0\",\n";
        //tFile << "{\n\t\"file-series-version\" : \"1.0\",\n";
        tString += "\t\"files\" : [\n";
        //tFile << "\t\"files\" : [\n";
        counter = 0;
    }
    //----------------------------------------------------
    void TSeriesWriteTime(std::string currFileName, prec time)
    {
        if (counter != 0)
        {
            tString += ",\n";
        }
        tString += "\t\t{ \"name\" : \"" + currFileName + "\", \"time\" : " + std::to_string(time) + "}";
        std::string temp_tString = tString + close_tString;
        if (binWrite)
        {
            tFile.open(path_Tseries_file, std::ios::out | std::ios::binary | std::ios::trunc);
        }
        else
        {
            tFile.open(path_Tseries_file, std::ios::out | std::ios::trunc);
        }
        tFile << temp_tString;
        counter++;
        tFile.close();
        //tFile << "\t\t{ \"name\" : \"" << currFileName << "\", \"time\" : " << time << "}";
    }
};