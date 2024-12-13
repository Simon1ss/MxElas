%% Results
%% P4-P3 Hu-Zhang elements (order is 5 and 4)
% 
% result =
% 
%   4¡Á7 table
% 
%     meshsize    sigma_L2error    order_sigma    u_L2error     order_u    assembletime    solvetime
%     ________    _____________    ___________    __________    _______    ____________    _________
% 
%          1          0.17841             0         0.051764         0        6.3723        0.03257 
%        0.5         0.007575        4.5578        0.0041666     3.635        7.5897        0.36641 
%       0.25       0.00028157        4.7497       0.00028162     3.887        13.935         3.9959 
%      0.125       9.5307e-06        4.8847       1.7958e-05    3.9711        123.93         130.79 
%


clear
close
clc

%% set lame constants, get exact solution and rhs
mu = 1/2;
ld = 1;
[param, sigma_exact, u_exact, f_exact] = smoothdata3D(mu,ld);

%% different mesh sizes
N1_ = [1,2,4,8]';
N_test = length(N1_);
sigma_L2error = zeros(N_test,1);
u_L2error = zeros(N_test,1);
assembletime = zeros(N_test,1);
solvetime = zeros(N_test,1);

%% test rate
for iter = 1 : N_test
    N1 = N1_(iter);

    %% generate mesh and basic data structure
    [node, elem, misc] = unitcubemesh(N1);

    tic;
    %% DOF for stress space
    elem2dof_sigma = dofHZ3DP4(elem, misc);
    elem2dof_u = dofDG3D(elem,3);
    DoFtensor = symtensor3HZP4(elem, misc);

    %% generate the stiffness matrix
    C = compliance3(node, elem, elem2dof_sigma, DoFtensor, param);
    B = div3(node, elem, elem2dof_sigma, elem2dof_u, DoFtensor);
    [Ndof_u,Ndof_sigma] = size(B);
    A = [C,B';B,sparse(Ndof_u,Ndof_u)];

    %% right hand side
    F20 = bulkIntegral3(node, elem, elem2dof_u, f_exact);
    F = [zeros(Ndof_sigma,1);-F20];
    assembletime(iter) = toc;

    %% solve the equations
    tic;
    SOLU = A\F;
    solvetime(iter) = toc;
    sigma_h = SOLU(1:Ndof_sigma);
    u_h = SOLU(Ndof_sigma+1:end);

    %% compute error
    sigma_L2error(iter) = stressError3(node, elem, elem2dof_sigma, DoFtensor, sigma_exact, sigma_h);
    u_L2error(iter) = uError3(node, elem, u_exact, u_h);
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

