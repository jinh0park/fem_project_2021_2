clear all;
close all;
clc;
L = 10;

%% Part 1. H

d_maxss = cell(1,3);
max_d_idxss = cell(1,3);
vonmises_maxss = cell(1,3);
max_vonmises_idxss = cell(1,3);
tresca_maxss = cell(1,3);
max_tresca_idxss = cell(1,3);
locationss = cell(1,3);
location_vonss = cell(1,3);
location_tress = cell(1,3);


x_hangers = [7 6 8];
ws = [2000 2330 1750];

for j=1:3
    w = ws(j);
    x_hanger = x_hangers(j);
    for i=1:20
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

    
        d_maxss{j} = d_maxs;
        max_d_idxss{j} = max_d_idxs;
        vonmises_maxss{j} = vonmises_maxs;
        max_vonmises_idxss{j} = max_vonmises_idxs;
        tresca_maxss{j} = tresca_maxs;
        max_tresca_idxss{j} = max_tresca_idxs;
        locationss{j} = locations;
        location_vonss{j} = location_vons;
        location_tress{j} = location_tres;
    
        fprintf("i=%d\n",i);
    end
    fprintf("j=%d\n",j);
end