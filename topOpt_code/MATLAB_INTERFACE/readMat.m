function [mat] = readMat(name)
     % open the file
     fid = fopen(char(name));
     if fid>0
        line = fgetl(fid);
        lineSet = textscan(line, "%d");
        nRow = lineSet{1}(1); nCol = lineSet{1}(2);
        mat = zeros(nRow, nCol);
        for i = 1:nRow
            line = fgetl(fid);
            lineSet = textscan(line, "%f");
            for j = 1:nCol
                mat(i,j) = lineSet{1}(j);
            end
        end
    end
    fclose(fid);
end