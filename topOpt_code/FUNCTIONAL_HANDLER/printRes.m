clc; 
clear; 
close all;

problemName = "double_pipe_extremely_fine";
test = "vel[1]_S_Vr[0.34]_goc_a1e5_b[1_1_0_0]_time[stat]_p2";
name = problemName + "/" + test;
% nameF = "optimization/" + name + "/func.txt";
% nameC = "optimization/" + name + "/changes.txt";
% nameV = "optimization/" + name + "/valid.txt";
nameF = name + "/func.txt";
name_no_w_F = name + "/no_weight_func.txt";
nameC = name + "/changes.txt";
nameV = name + "/valid.txt";
nfig = 0;

print_func = 1;
print_changes = 1;
print_no_weight_func = 1;

if (print_func == 1)
    functional = readMat(nameF);
    f0Init = functional(1,1);
    f0Init_out_box = functional(1,2);
    f0Init_parts = functional(1,3:end);
    functional = functional(2:end, :);
    func = functional(1:end,1);
    func_out_box = functional(1:end,2);
    func_alpha = functional(1:end,3);
    func_grad  = functional(1:end,4);
    func_vort  = functional(1:end,5);
    func_p     = functional(1:end,6);
    
    valid = discardVec(nameV);
    
    % fid = fopen("optimization/" + name + "/latex.txt", 'wt');
    % fid = fopen(name + "/latex.txt", 'wt');
    % 
    % for i = 1:length(func)
    %     fprintf(fid, "(%d, %f) \n", i, func(i));
    % end
    % fclose(fid);
    
    it = 1:length(func);
    nfig = nfig + 1;
    figure(nfig);
    hold on;
    plot(it, func, "-", "Color", "#000000", "LineWidth", 2);
    plot(it, func_out_box,"-", "Color", "#A2142F", "LineWidth", 2);
    plot(it, func_alpha,"-", "Color", "#D95319", "LineWidth", 2);
    plot(it, func_grad,"-", "Color", "#0072BD", "LineWidth", 2);
    plot(it, func_vort,"-", "Color", "#77AC30", "LineWidth", 2);
    plot(it, func_p, "-", "Color", "#EDB120", "LineWidth", 2);
    for i = 1 : floor(length(func) / 20) : (lengtchangesh(func))
        if valid(i) == 1
            plot(it(i), func(i), 'bo', "LineWidth", 1.5);
        else
            plot(it(i), func(i), 'ro', "LineWidth", 1.5);
        end
    end
    legend("Jtot", "Jout", "Jalpha", "Jgrad", "Jvort", "Jp");
    title("FUNCTIONAL");
    % xscale("log");
    hold off;
end
 

if (print_changes == 1)    
    changes =  discardVec(nameC);
    nfig = nfig + 1;
    figure(nfig);
    
    plot(it, changes, "o-");
    
    title("CHANGES");
end

if (print_no_weight_func == 1)    
    no_w_functional = readMat(name_no_w_F);
    no_w_f_alpha = no_w_functional(:,1);
    no_w_f_grad = no_w_functional(:,2);
    no_w_f_vort = no_w_functional(:,3);
    no_w_f_p = no_w_functional(:,4);

    nfig = nfig + 1;
    figure(nfig);
    hold on
    plot(it, no_w_f_alpha,"-", "Color", "#D95319", "LineWidth", 2);
    plot(it, no_w_f_grad,"-", "Color", "#0072BD", "LineWidth", 2);
    plot(it, no_w_f_vort,"-", "Color", "#77AC30", "LineWidth", 2);
    plot(it, no_w_f_p, "-", "Color", "#EDB120", "LineWidth", 2);
    legend("Jalpha", "Jgrad", "Jvort", "Jp");
    title("NO W FUNCTIONAL");
    hold off;
end




