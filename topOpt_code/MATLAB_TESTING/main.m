clc; clear; close all;

mat00 = discard("mat00.txt");
mat01 = discard("mat01.txt");
mat02 = discard("mat02.txt");
mat03 = discard("mat03.txt");
mat10 = discard("mat10.txt");
mat11 = discard("mat11.txt");
mat12 = discard("mat12.txt");
mat13 = discard("mat13.txt");
mat20 = discard("mat20.txt");
mat21 = discard("mat21.txt");
mat22 = discard("mat22.txt");
mat23 = discard("mat23.txt");

mat = discard("mat.txt");

matMat = [mat00, mat01, mat02, mat03; mat10, mat11, mat12, mat13; mat20, mat21, mat22, mat23];

diff = matMat - mat;

norm(diff, "fro")