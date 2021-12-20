clear all;
close all;
clc;
L = 10;
%% Part 1.E
x_hanger = 7;
w = 2000;
M = 10;


[d, sigmas] = fem_triangle_function(x_hanger, w, M);

[d_max, max_idx] = max(abs(d));

fprintf("Node %d has maximum displacement: %d [m].\n which is " + ...
    "%.2f m far from left end.\n", ...
    ceil(max_idx/2), d_max, ceil((max_idx-1)/4)*L/M);

d_y_upper = d(2:4:end);
d_y_lower = d(4:4:end);

plot(d_y_lower*10^7);
hold on
plot(d_y_upper*10^7+1);
hold off
%% Part 1.F
pr_sigmas = zeros(3, 2*M);

pr_sigmas(1,:) = (sigmas(1,:) + sigmas(2,:))/2 + sqrt(((sigmas(1,:)-sigmas(2,:))/2).^2+sigmas(3,:).^2);
pr_sigmas(2,:) = (sigmas(1,:) + sigmas(2,:))/2 - sqrt(((sigmas(1,:)-sigmas(2,:))/2).^2+sigmas(3,:).^2);

von_mises = sqrt(((pr_sigmas(1,:)-pr_sigmas(2,:)).^2+(pr_sigmas(2,:)-pr_sigmas(3,:)).^2+(pr_sigmas(3,:)-pr_sigmas(1,:)).^2)/2);
tresca = max(pr_sigmas) - min(pr_sigmas);