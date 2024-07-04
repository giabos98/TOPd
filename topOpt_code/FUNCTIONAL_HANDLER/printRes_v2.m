clc; 
clear; 
close all;

nfig = 1;

comparison = Compare_Functs();

problemName = "dp_square_inv_1_reg_d5";
test1 = "case1_Re[1]_B[100]";
test2 = "case1_Re[1]_B[1000]";
test3 = "case1_Re[1]_B[10000]";
test4 = "case1_Re[1]_B[10000]_smooth";
% test5 = "rho[1]_mu[1]_U[10]_a[10e4]_L[1]_mma";
% test6 = "rho[1]_mu[1]_U[10]_a[10e4]_L[1]_goc";
% % test7 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_mma";
% % test8 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_goc";

functional1 = Functional(problemName, test1, 1, "beta 100   ");
functional2 = Functional(problemName, test2, 1, "beta 1000  ");
functional3 = Functional(problemName, test3, 1, "beta 10000 ");
functional4 = Functional(problemName, test4, 1, "beta smooth");
% functional5 = Functional(problemName, test5, 1, " 10 MMA");
% functional6 = Functional(problemName, test6, 1, " 10 GOC");
% % functional7 = Functional(problemName, test7, 1, "100 MMA");
% % functional8 = Functional(problemName, test8, 1, "100 GOC");

comparison = comparison.add_functional(functional1);
comparison = comparison.add_functional(functional2);
comparison = comparison.add_functional(functional3);
comparison = comparison.add_functional(functional4);
% comparison = comparison.add_functional(functional5);
% comparison = comparison.add_functional(functional6);
% comparison = comparison.add_functional(functional7);
% comparison = comparison.add_functional(functional8);


% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
comparison.compare_functionals(nfig, "grad", 0, 0, "", [], 300);
% comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);





