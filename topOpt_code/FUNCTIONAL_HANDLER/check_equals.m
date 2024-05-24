clc; clear; close all;

test_file_name1 = "test1";
path1 = "../scratch/" + test_file_name1 + ".txt";

test_file_name2 = "test2";
path2 = "../scratch/" + test_file_name2 + ".txt";

mat1 = readMat(path1);
mat2 = readMat(path2);
test_len = size(mat1,1);

for i = 1:test_len
    if (mat1(i,1) ~= mat2(i,1))
        fprintf("comp1: error in " + int2str(i));
    elif (mat1(i,2) ~= mat2(i,2))
        fprintf("comp2: error in " + int2str(i));
    else
        % fprintf("TRUE\n");
    end
end