%% Results
%% P3-P2 Hu-Zhang elements (order is 4 and 3)
% result =
% 
%   6¡Á7 table
% 
%     meshsize    sigma_L2error    order_sigma    u_L2error     order_u    assembletime    solvetime
%     ________    _____________    ___________    __________    _______    ____________    _________
% 
%          0.5      0.0029597             0       0.00088976         0       0.043247       0.001105
%         0.25     0.00026152        3.5005       0.00013926    2.6756       0.048611      0.0057867
%        0.125      1.878e-05        3.7997        1.848e-05    2.9138       0.065982       0.029793
%       0.0625     1.2434e-06        3.9168       2.3462e-06    2.9776        0.13606        0.18196
%      0.03125     7.9719e-08        3.9632       2.9443e-07    2.9943        0.49586         1.1974
%     0.015625     5.0415e-09         3.983        3.684e-08    2.9986         2.6112         11.434
%


clear
close all
clc

%% set lame constants, get exact solution and rhs
mu = 1/2;
ld = 1;
[param, sigma_exact, u_exact, f_exact] = smoothdata2D(mu,ld);

%% different mesh sizes
N1_ = [2,4,8,16,32,64]';
N_test = length(N1_);
sigma_L2error = zeros(N_test,1);
u_L2error = zeros(N_test,1);
assembletime = zeros(N_test,1);
solvetime = zeros(N_test,1);

%% test rate
for iter = 1 : N_test
    N1 = N1_(iter);
    
    %% generate mesh and basic data structure
    [node, elem, ~] = unitsquaremesh(N1);
    
    tic;
    %% DOF
    elem2dof_sigma = dofHZ2DP3(elem);
    elem2dof_u = dofDG2D(elem,2);
    DoFtensor = symtensor2HZP3(elem,node);
    
    %% matrices
    C = compliance2(node,elem,elem2dof_sigma,DoFtensor,param);
    B = div2(node,elem,elem2dof_sigma,elem2dof_u,DoFtensor);
    [Ndof_u,Ndof_sigma] = size(B);
    A = [C,B';B,sparse(Ndof_u,Ndof_u)];
    
    %% right hand side
    F20 = bulkIntegral2(node, elem, elem2dof_u, f_exact);
    F = [zeros(Ndof_sigma,1);-F20];
    assembletime(iter) = toc;

    %% solve the equations
    tic;
    SOLU = A\F;
    solvetime(iter) = toc;
    sigma_h = SOLU(1:Ndof_sigma);
    u_h = SOLU(Ndof_sigma+1:end);
    
    %% compute L2 error
    sigma_L2error(iter) = stressError2(node, elem, elem2dof_sigma, DoFtensor, sigma_exact, sigma_h);
    u_L2error(iter) = uError2(node, elem, u_exact, u_h);
end



%% display convergence order
meshsize = 1./N1_;
order_sigma = zeros(N_test,1);
order_u = zeros(N_test,1);
for iter = 2 : N_test
    order_sigma(iter) = log(sigma_L2error(iter-1)/sigma_L2error(iter))/...
        log(N1_(iter)/N1_(iter-1));
    order_u(iter) = log(u_L2error(iter-1)/u_L2error(iter))/...
        log(N1_(iter)/N1_(iter-1));
end
result = table(meshsize,sigma_L2error,order_sigma,u_L2error,order_u,assembletime,solvetime);
display(result);


figure(1)
showrateh(meshsize,sigma_L2error);
figure(2)
showrateh(meshsize,u_L2error);
