clc; 
clear; 
close all;

nfig = 1;

comparison_square = Compare_Functs();
comparison_rect = Compare_Functs();

problemName = "bifuraction_L1.5_d10";
test1 = "case1_Re[1]_B[1e4]_smooth";
test2 = "case1_Re[1]_B[1e4]_smooth_convex";
test3 = "case1_Re[1]_B[1e4]_smooth_concave";
% test4 = "rect_Re[1]_B[100]";
% test5 = "rect_Re[1]_B[10000]";
% test6 = "rect_Re[1]_B[10000]_smooth";
% % test7 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_mma";
% % test8 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_goc";

functional1 = Functional(problemName, test1, 1, "linear");
functional2 = Functional(problemName, test2, 1, "convex");
functional3 = Functional(problemName, test3, 1, "concave");
% functional4 = Functional(problemName, test4, 1, "beta 1e2   ");
% functional5 = Functional(problemName, test5, 1, "beta 1e4  ");
% functional6 = Functional(problemName, test6, 1, "beta 1e4 (smooth) ");
% % functional7 = Functional(problemName, test7, 1, "100 MMA");
% % functional8 = Functional(problemName, test8, 1, "100 GOC");

comparison_square = comparison_square.add_functional(functional1);
comparison_square = comparison_square.add_functional(functional2);
comparison_square = comparison_square.add_functional(functional3);
% comparison_rect = comparison_rect.add_functional(functional4);
% comparison_rect = comparison_rect.add_functional(functional5);
% comparison_rect = comparison_rect.add_functional(functional6);
% comparison = comparison.add_functional(functional7);
% comparison = comparison.add_functional(functional8);


colors = ["#00B050", "#C00000", "#00B0F0"];

% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
nfig = comparison_square.compare_functionals(nfig, "tot", 0, 1, "", colors, 200);
% comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);7



