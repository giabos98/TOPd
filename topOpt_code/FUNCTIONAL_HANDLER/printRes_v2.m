clc; 
clear; 
close all;

nfig = 1;

comparison = Compare_Functs();

problemName = "double_pipe_reg_0.1";
test1 = "rho[1]_mu[1]_U[0.1]_a[0.1e4]_L[1]_mma";
test2 = "rho[1]_mu[1]_U[0.1]_a[0.1e4]_L[1]_goc";
test3 = "rho[1]_mu[1]_U[1]_a[1e4]_L[1]_mma";
test4 = "rho[1]_mu[1]_U[1]_a[1e4]_L[1]_goc";
test5 = "rho[1]_mu[1]_U[10]_a[10e4]_L[1]_mma";
test6 = "rho[1]_mu[1]_U[10]_a[10e4]_L[1]_goc";
% test7 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_mma";
% test8 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_goc";

functional1 = Functional(problemName, test1, 1, "0.1 MMA");
functional2 = Functional(problemName, test2, 1, "0.1 GOC");
functional3 = Functional(problemName, test3, 1, "  1 MMA");
functional4 = Functional(problemName, test4, 1, "  1 GOC");
functional5 = Functional(problemName, test5, 1, " 10 MMA");
functional6 = Functional(problemName, test6, 1, " 10 GOC");
% functional7 = Functional(problemName, test7, 1, "100 MMA");
% functional8 = Functional(problemName, test8, 1, "100 GOC");

comparison = comparison.add_functional(functional1);
comparison = comparison.add_functional(functional2);
comparison = comparison.add_functional(functional3);
comparison = comparison.add_functional(functional4);
comparison = comparison.add_functional(functional5);
comparison = comparison.add_functional(functional6);
% comparison = comparison.add_functional(functional7);
% comparison = comparison.add_functional(functional8);


% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
% comparison.compare_functionals(nfig, "tot", 0, 0, "", []);
comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);





