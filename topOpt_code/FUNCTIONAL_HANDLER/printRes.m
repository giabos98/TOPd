clc; 
clear; 
close all;

nfig = 1;

comparison_square = Compare_Functs();
comparison_rect = Compare_Functs();

problemName = "dp_inv_1_reg_d4";
test1 = "test_power_law";
% test2 = "rect_Re[1]_B[10000]";
% test3 = "rect_Re[1]_B[10000]_smooth";
% test4 = "rect_Re[1]_B[1e4]_smooth_line";
% test5 = "rect_Re[1]_B[1e4]_smooth_convex";
% test6 = "rect_Re[1]_B[1e4]_smooth_concave";
% % test7 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_mma";
% % test8 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_goc";

functional1 = Functional(problemName, test1, 1, "");
% functional2 = Functional(problemName, test2, 1, "1e4");
% functional3 = Functional(problemName, test3, 1, "1e4 smooth");
% functional4 = Functional(problemName, test4, 1, "1e4 linear");
% functional5 = Functional(problemName, test5, 1, "1e4 convex");
% functional6 = Functional(problemName, test6, 1, "1e4 concave");
% % functional7 = Functional(problemName, test7, 1, "100 MMA");
% % functional8 = Functional(problemName, test8, 1, "100 GOC");

comparison_rect = comparison_rect.add_functional(functional1);
% comparison_rect = comparison_rect.add_functional(functional2);
% comparison_rect = comparison_rect.add_functional(functional3);
% comparison_rect = comparison_rect.add_functional(functional4);
% comparison_rect = comparison_rect.add_functional(functional5);
% comparison_rect = comparison_rect.add_functional(functional6);
% comparison = comparison.add_functional(functional7);
% comparison = comparison.add_functional(functional8);


% colors = ["#00B050", "#C00000", "#00B0F0"];
colors = [];

% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
nfig = comparison_rect.compare_functionals(nfig, "tot", 0, 1, "", colors, 200);
nfig = comparison_rect.compare_alpha_max(nfig);
% comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);7



