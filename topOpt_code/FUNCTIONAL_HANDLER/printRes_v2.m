clc; 
clear; 
close all;

nfig = 1;

comparison = Compare_Functs();

problemName = "double_pipe_reg_0.02";
% test1 = "vel[0.9]_[S]_Re[0.1]_Vr[0.34]_[goc]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
% test2 = "vel[0.9]_[S]_Re[0.1]_Vr[0.34]_[mma]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test3 = "vel[0.9]_[NS]_Re[0.1]_Vr[0.34]_[goc]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test4 = "vel[0.9]_[NS]_Re[0.1]_Vr[0.34]_[mma]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test5 = "vel[0.9]_[NS]_Re[1]_Vr[0.34]_[goc]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test6 = "vel[0.9]_[NS]_Re[1]_Vr[0.34]_[mma]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test7 = "vel[0.9]_[NS]_Re[10]_Vr[0.34]_[goc]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
test8 = "vel[0.9]_[NS]_Re[10]_Vr[0.34]_[mma]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
% test9 = "vel[0.9]_[NS]_Re[100]_Vr[0.34]_[goc]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";
% test10 = "vel[0.9]_[NS]_Re[100]_Vr[0.34]_[mma]_a[1e4]_b[1_1_0_0]_t[stat]_p[2]_m[1]_s[0]";

% functional1 = Functional(problemName, test1, "S GOC");
% functional2 = Functional(problemName, test2, "S MMA");
functional3 = Functional(problemName, test3, 1, "NS Re=1 GOC");
functional4 = Functional(problemName, test4, 1, "NS Re=0.1 MMA");
functional5 = Functional(problemName, test5, 0.3,"NS Re=1 GOC");
functional6 = Functional(problemName, test6, 0.3,"NS Re=1 MMA");
functional7 = Functional(problemName, test7, 0.1,"NS Re=10 GOC");
functional8 = Functional(problemName, test8, 0.1, "NS Re=10 MMA");
% functional9 = Functional(problemName, test9, 1, "NS Re=100 GOC");
% functional10 = Functional(problemName, test10, 1, "NS Re=100 MMA");

% comparison = comparison.add_functional(functional1);
% comparison = comparison.add_functional(functional2);
comparison = comparison.add_functional(functional3);
comparison = comparison.add_functional(functional4);
comparison = comparison.add_functional(functional5);
comparison = comparison.add_functional(functional6);
comparison = comparison.add_functional(functional7);
comparison = comparison.add_functional(functional8);
% comparison = comparison.add_functional(functional9);
% comparison = comparison.add_functional(functional10);


% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
comparison.compare_functionals(nfig, "tot", 0, "", []);
% nfig = comparison.print(nfig, [1,0,0]);





