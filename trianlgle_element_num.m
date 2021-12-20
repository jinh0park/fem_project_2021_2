clear all;
close all;
clc;
L = 10;
%% Part 1. G
x_hanger = 7;
w = 2000;

for i=20:20
    M = i*10; % Element Number  = k = 2*M = [20, 40, 60, ... , 400]

    [d, sigmas] = fem_triangle_function(x_hanger, w, M);
    [d_max, max_d_idx] = max(abs(d));

    % Stress
    pr_sigmas = zeros(3, 2*M);
    
    pr_sigmas(1,:) = (sigmas(1,:) + sigmas(2,:))/2 + sqrt(((sigmas(1,:)-sigmas(2,:))/2).^2+sigmas(3,:).^2);
    pr_sigmas(2,:) = (sigmas(1,:) + sigmas(2,:))/2 - sqrt(((sigmas(1,:)-sigmas(2,:))/2).^2+sigmas(3,:).^2);
    
    von_mises = sqrt(((pr_sigmas(1,:)-pr_sigmas(2,:)).^2+(pr_sigmas(2,:)-pr_sigmas(3,:)).^2+(pr_sigmas(3,:)-pr_sigmas(1,:)).^2)/2);
    tresca = max(pr_sigmas) - min(pr_sigmas);

    [vonmises_max, max_vonmises_idx] = max(von_mises);
    [tresca_max, max_tresca_idx] = max(tresca);

    d_maxs(i) = d_max;
    max_d_idxs(i) = max_d_idx;
    vonmises_maxs(i) = vonmises_max;
    max_vonmises_idxs(i) = max_vonmises_idx;
    tresca_maxs(i) = tresca_max;
    max_tresca_idxs(i) = max_tresca_idx;
    locations(i) = floor((max_d_idx-1)/4)*L/M;
    location_vons(i) = floor((max_vonmises_idx-1)/2)*L/M + L/(2*M);
    location_tres(i) = floor((max_tresca_idx-1)/2)*L/M + L/(2*M);

    fprintf("i=%d\n",i);
end