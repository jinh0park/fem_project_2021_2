%% 2017-15024 기계공학부 박진호
clear all;
close all;
clc;
%% Constants
L=10;           % Length of the girder [m]
b=1;            % Width of the girder [m]
t=0.5;          % Thickness of the girder [m]
x_hanger=7;     % Length of the region which is loaded by hanger [m]
w=2*10^3;       % Distributed load due to hanger [N/m^2]
v=0.3;          % Poisson ration [1]
E=200*10^9;     % Youngs modulus of the girder [Pa]
Y=800*10^6;     % Yield Strength of the girder [Pa]
%% Define Varialbes
M=10;              % The number of divided section of the length (=k/2)
k=M;           % The number of rectangle elements
n=2*k+2;          % The number of nodes while using triangle elements
D=E/(1-v^2)*[1 v 0; 
                 v 1 0; 
                 0 0 (1-v)/2];   % D matrix for plane stress
%% Set (x, y) coordinates of the nodes
% 1 - 3 - 5 - ... -(n-1)  
% |   |   |   ...   |
% 2 - 4 - 6 - ... - n
x_i=zeros(1,n);
y_i=zeros(1,n);
for i=1:M+1
    x_i(2*i-1)=L*(i-1)/M;
    x_i(2*i)=L*(i-1)/M;
    y_i(2*i-1)=b;
    y_i(2*i)=0;
end
%% Find the node numbers of the each rectangle elements
element=zeros(4, k);
for i=1:k
    ii=2*i-1;
    element(:,i)=[ii ii+1 ii+3 ii+2];
end
%% Boundary Conditions
u_first=ceil(M*(L-x_hanger)/(2*L))+1;   % First element location which is loaded
u_last=M+2-u_first;                     % Last element location which is loaded  
f=(w*t*x_hanger)/(u_last-u_first);
F=zeros(1,2*n);
for u=u_first:u_last
    F(4*u)=-f;
end
F_known=F(5:2*n-4); % F without each ends of the girder
                    % which is known as boundary conditions.
%% Construct Connectivity vectors
% Neccessay to add local K matrices to global K matrix.
C=cell(k,1);
for i=1:k
    c=zeros(8, 2*n);
    for j=1:4
        c(2*j-1:2*j, 2*element(j,i)-1:2*element(j,i))=eye(2);
    end
    C{i}=c;
end
%% Define Shape Functions
syms x y r s;       % For isoparametric method, (x, y) -> (r, s)
N(1)=(1-r)*(1+s)/4;
N(2)=(1-r)*(1-s)/4;
N(3)=(1+r)*(1-s)/4;
N(4)=(1+r)*(1+s)/4;
%% Set Global Stiffness Matrix K
K=zeros(2*n);       % DOF = 2 * node number(=n)
%% Build Local Stiffness Matrix k
Bs=cell(k,1);       % Store B to calculate stress later.
for element_index=1:k
    i=element(1, element_index);
    j=element(2, element_index);
    m=element(3, element_index);
    q=element(4, element_index);
    % Compute Jacobian matrix
    x=x_i(i)*N(1)+x_i(j)*N(2)+x_i(m)*N(3)+x_i(q)*N(4);
    y=y_i(i)*N(1)+y_i(j)*N(2)+y_i(m)*N(3)+y_i(q)*N(4);
    J=jacobian([x y], [r s]);
    % Compute B matrix and local
    B=sym(zeros(3,8));
    for p=1:4
        B_i=zeros(3,2);
        dN=inv(J)*[diff(N(p),r) diff(N(p),s)].';
        B_i=[dN(1) 0; 0 dN(2); dN(2) dN(1)];
        B(:,2*p-1:2*p)=B_i;
    end
    Bs{element_index}=B;
    % Compute local k matrix
    K_local=double(int(int(B'*D*B*t*det(J),r,-1,1),s,-1,1));
    % Add to Global K matrix
    K=K+C{element_index}'*K_local*C{element_index};
end
%% Calculate Displacements and Forces
K_part=K(5:2*n-4, 5:2*n-4);
d=zeros(1,2*n);
d(5:2*n-4)=inv(K_part)*F_known';
F=K*d';
%% Answer - Maximum Displacement
[d_max, idx]=max(abs(d));
fprintf("Node %d has maximum displacement: %d [m].\n which is " + ...
    "%.2f m far from left end.\n", ...
    ceil(idx/2), d_max, ceil((idx-1)/4)*L/M);
%% Answer - Maximum Stress
sigmas=zeros(3,k);
for element_index=1:k
    i=element(1, element_index);
    j=element(2, element_index);
    m=element(3, element_index);
    q=element(4, element_index);
    d_partial=zeros(1,8);
    d_partial(1:2)=d(2*element(1,element_index)-1:2*element(1,element_index));
    d_partial(3:4)=d(2*element(2,element_index)-1:2*element(2,element_index));
    d_partial(5:6)=d(2*element(3,element_index)-1:2*element(3,element_index));
    d_partial(7:8)=d(2*element(4,element_index)-1:2*element(4,element_index));
    sigma = D*Bs{element_index}*d_partial'; % Stress of an element was found
    sigmas(:,element_index)=subs(sigma, [r, s], [0, 0]);
end
