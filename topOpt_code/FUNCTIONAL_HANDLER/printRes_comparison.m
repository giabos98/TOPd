clc; 
clear; 
close all;

nfig = 1;

comparison_square = Compare_Functs();
comparison_rect = Compare_Functs();

problemName = "invariance_results";
test1 = "square_Re[1]_B[100]";
test2 = "square_Re[1]_B[10000]";
test3 = "square_Re[1]_B[10000]_smooth";
test4 = "rect_Re[1]_B[100]";
test5 = "rect_Re[1]_B[10000]";
test6 = "rect_Re[1]_B[10000]_smooth";
% % test7 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_mma";
% % test8 = "rho[1]_mu[1]_U[100]_a[100e4]_L[1]_goc";

functional1 = Functional(problemName, test1, 1, "beta 1e2   ");
functional2 = Functional(problemName, test2, 1, "beta 1e4  ");
functional3 = Functional(problemName, test3, 1, "beta 1e4 (smooth) ");
functional4 = Functional(problemName, test4, 1, "beta 1e2   ");
functional5 = Functional(problemName, test5, 1, "beta 1e4  ");
functional6 = Functional(problemName, test6, 1, "beta 1e4 (smooth) ");
% % functional7 = Functional(problemName, test7, 1, "100 MMA");
% % functional8 = Functional(problemName, test8, 1, "100 GOC");

comparison_square = comparison_square.add_functional(functional1);
comparison_square = comparison_square.add_functional(functional2);
comparison_square = comparison_square.add_functional(functional3);
comparison_rect = comparison_rect.add_functional(functional4);
comparison_rect = comparison_rect.add_functional(functional5);
comparison_rect = comparison_rect.add_functional(functional6);
% comparison = comparison.add_functional(functional7);
% comparison = comparison.add_functional(functional8);


colors = ["#00B050", "#C00000", "#00B0F0"];

% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
nfig = comparison_square.compare_functionals(nfig, "tot", 0, 1, "", colors, 200);
% comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);7

figure(nfig);
hold on;
xscale("log");
yscale("log");
xticks([1 10 50 100 200]);
yticks([1:2:11]);
axis([0.9 200 0.9 11]);
set(gca,'linewidth',2);
set(gca,'fontsize', 28, "FontName", "Times");
% xlabel("Iterations", "FontSize", 40, "FontName", "Times");
% ylabel("Relative Functional", "FontSize", 36);
hold off;


% comparison.compare_base_functionals(nfig, [-2,1], 0, "", "tot");
nfig = comparison_rect.compare_functionals(nfig, "tot", 0, 1, "", colors, 200);
% comparison.compare_functionals_initial_values(nfig);
% nfig = comparison.print(nfig, [1,0,0]);
figure(nfig);
hold on;
xscale("log");
yscale("log");
xticks([1 10 50 100 200]);
yticks([1,2:2:14]);
axis([0.9 200 0.9 14]);
set(gca,'linewidth',2);
set(gca,'fontsize', 28, "FontName", "Times");
% xlabel("Iterations", "FontSize", 40, "FontName", "Times");
% ylabel("Relative Functional", "FontSize", 36);
hold off;


%% PRINT ALPHA_MAX(IT)
nfig = nfig+1;
% it = [[1, 10]; [10, 20]; [20, 100]];
% alpha_max = [[1, 1e2]; [1e2,1e3]; [1e3, 1e4]];
n = 4;
figure(nfig);
hold on;
it = [];
beta = [];
it = 1:0.1:100;
beta_line(1:991) = 1e2 + (1e4-1e2)*(it(1:end)-1)/(it(end)-1);
beta_convex(1:491) = 1e2 + (1e3-1e2)*(it(1:491)-1)/(50-1);
beta_convex(492:991) = 1e3 + (1e4-1e3)*(it(492:end)-50)/(100-50);
beta_concave(1:491) = 1e2 + (9.1e3-1e2)*(it(1:491)-1)/(50-1);
beta_concave(492:991) = 9.1e3 + (1e4-9.1e3)*(it(492:end)-50)/(100-50);
% plot(it, beta,"-k", "LineWidth", 5);
z = zeros(size(it));
green = hex2rgb(colors(1));
red = hex2rgb(colors(2));
yellow = hex2rgb("#FFC000");
map = [];
% for i=1:1:81
%     map = [map; green];
% end

delta_it = 491;
for i=1:1:delta_it
    rel_it = (i-1)/delta_it;
    col = green*(1-rel_it) + yellow*rel_it;
    map = [map; col];
end

delta_it = length(it)-(delta_it+1);
for i=delta_it+1:1:length(it)
    rel_it = (i-(delta_it+1))/delta_it;
    col = yellow*(1-rel_it) + red*rel_it;
    map = [map; col];
end
small_dot_line = 4;
big_dot_line  = 5;
small_marker_size = 28;
big_marker_size = 34;
big_line_width = 24;
small_line_width = 18;
surf([it;it],[beta_line;beta_line],[z;z],[it;it],'facecol','no','edgecol','interp','linew', big_line_width);
surf([it;it],[beta_convex;beta_convex],[z;z],[it;it],'facecol','no','edgecol','interp','linew', big_line_width);
surf([it;it],[beta_concave;beta_concave],[z;z],[it;it],'facecol','no','edgecol','interp','linew', big_line_width);
% plot(1:0.1:9, ones(81,1)*1e2, "-", "Color", green, "LineWidth", small_line_width);
plot(101:1:110, ones(10,1)*1e4, "-", "Color", red, "LineWidth", small_line_width);
plot(111:1:130, ones(20,1)*1e4, ":", "Color", red, "LineWidth", small_line_width);
plot(1, 1e2, 'ko', "LineWidth", small_dot_line, "MarkerSize", small_marker_size, 'MarkerFaceColor',green);
% plot(10, 1e2, 'ko', "LineWidth", big_dot_line, "MarkerSize", big_marker_size, 'MarkerFaceColor',green);
k_50 = find(it==50);
plot(50, beta_line(k_50), 'ko', "LineWidth", big_dot_line, "MarkerSize", big_marker_size, 'MarkerFaceColor', yellow);
plot(50, beta_convex(k_50), 'ko', "LineWidth", big_dot_line, "MarkerSize", big_marker_size, 'MarkerFaceColor', yellow);
plot(50, beta_concave(k_50), 'ko', "LineWidth", big_dot_line, "MarkerSize", big_marker_size, 'MarkerFaceColor', yellow);
% k = find(it==40);
% plot(40, beta(k), 'ko', "LineWidth", small_dot_line, "MarkerSize", small_marker_size, 'MarkerFaceColor', map(k,:));
% k = find(it==60);
% plot(60, beta(k), 'ko', "LineWidth", small_dot_line, "MarkerSize", small_marker_size, 'MarkerFaceColor', map(k,:));
% k = find(it==80);
% plot(80, beta(k), 'ko', "LineWidth", small_dot_line, "MarkerSize", small_marker_size, 'MarkerFaceColor', map(k,:));
plot(100, 1e4, 'ko', "LineWidth", big_dot_line, "MarkerSize", big_marker_size, 'MarkerFaceColor', red);
% plot(115, 1e4, 'ko', "LineWidth", small_dot_line, "MarkerSize", small_marker_size, 'MarkerFaceColor', red);

colormap(map);
% colorbar;
% xscale("log");
% yscale("log");
xticks([1 10 50 100]);
yticks([1e2 beta_convex(k_50) beta_line(k_50) beta_concave(k_50) 1e4]);
yticklabels(["10^2" "10^3" "5x10^3" "9.1x10^3" "10^4"]);
axis([-5 130 -600 1e4]);
set(gca,'linewidth',small_dot_line);
set(gca,'fontsize', 48, "FontName", "Times");
set(gca,'Color','None')
% grid on;
% xlabel("Iterations \theta", "FontSize", 40, "FontName", "Times");
% ylabel("\beta(\theta)", "FontSize", 36);
hold off;


