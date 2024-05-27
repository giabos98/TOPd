clc; 
clear; 
close all;

nfig = 1;

comparison = Compare_Functs();

problemName = "double_pipe_extremely_fine";
test1 = "vel[1]_S_Vr[0.34]_mma_a1e5_b[1_1_0_0]_time[stat]_p2";
test2 = "vel[1]_S_Vr[0.34]_goc_a1e5_b[1_1_0_0]_time[stat]_p2";
test3 = "vel[1]_S_Vr[0.34]_mma_a1e5_b[1_1_0_0]_time[stat]_p2";
test4 = "vel[1]_S_Vr[0.34]_goc_a1e5_b[1_1_0_0]_time[stat]_p2";



functional1 = Functional(problemName, test1, "mma_S");
functional2 = Functional(problemName, test2, "goc_S");

comparison = comparison.add_functional(functional1);
comparison = comparison.add_functional(functional2);

comparison.compare_base_functionals(nfig);
% nfig = comparison.print(nfig, [1,0,0]);





