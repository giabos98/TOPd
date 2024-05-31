clc; close all; clear;

nfig = 0;

l = 0.1:0.001:10;

mu = 0.001:0.01:10;

s = 0;
for i = 1:length(l)
    if (l(i)==1)
        s = i;
    end
end

rho = zeros(length(mu), length(l));
alpha = zeros(length(mu), length(l));

[X,Y]=meshgrid(mu,l);
x = X;
y = X./Y;
z = X./(Y.^2);


nfig = nfig+1;
figure(nfig);
sl = surf(x,y,z);
hold on;
plot3(x(1,:),y(1,:), z(1,:), "-k", "LineWidth", 2);
plot3(x(length(l),:),y(length(l),:), z(length(l),:), "-k", "LineWidth", 2);
plot3(x(:,length(mu)),y(:,length(mu)), z(:,length(mu)), "-k", "LineWidth", 2);
plot3(x(s,:),y(s,:), z(s,:), "-", "LineWidth", 4, "Color", "#EDB120");

colormap summer
sl.EdgeColor = 'none';
hold off;
