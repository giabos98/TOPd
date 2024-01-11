function [vec] = discardVec(name)
     % open the file
     fid = fopen(char(name));
     vec = zeros(100,1);     
     count = 0;
     if fid>0
         while ~feof(fid)
            count = count + 1;
            line = fgetl(fid);
            lineSet = textscan(line, "%f");
            vec(count) = lineSet{1}(1);
         end
    end
    fclose(fid);
    vec = vec(1:count);
end