function [mat] = discard(name)
     % open the file
     fid = fopen(char(name));
     
     i_sparse = zeros(80000,1);
     j_sparse = zeros(80000,1);
     coef = zeros(80000,1);
     
     count = 0;
     if fid>0
         while ~feof(fid)
            count = count + 1;
            line = fgetl(fid);
            lineSet = textscan(line, "%f");
            i_sparse(count) = lineSet{1}(1)+1;
            j_sparse(count) = lineSet{1}(2)+1;
            coef(count)     = lineSet{1}(3);
         end
    end
    fclose(fid);
    i_sparse = i_sparse(1:count); 
    j_sparse = j_sparse(1:count);
    coef = coef(1:count);
    mat = sparse(i_sparse, j_sparse, coef);
end