%% Results
%% P3-P2 Hu-Zhang elements (order is 4 and 3)
%% mixed boundary conditions
%   corner points between Dirichlet boundary and Neumann boundary might be "over-constrained"
% 
% result =
% 
%   6×7 table
% 
%     meshsize    sigma_L2error    order_sigma    u_L2error     order_u    assembletime    solvetime
%     ________    _____________    ___________    __________    _______    ____________    _________
% 
%          0.5      0.0075158             0          0.00658         0        0.10407      0.0015308
%         0.25     0.00057975        3.6964       0.00082856    2.9894       0.056924        0.00691
%        0.125     3.8195e-05         3.924       0.00010375    2.9975       0.070236       0.032585
%       0.0625     2.4027e-06        3.9906       1.2976e-05    2.9993        0.13108         0.1595
%      0.03125     1.4949e-07        4.0066       1.6222e-06    2.9998        0.43195        0.92069
%     0.015625     9.2978e-09         4.007       2.0278e-07         3          2.542          8.762


clear
close all
clc

%% set lame constants, get exact solution and rhs
mu = 1/2;
ld = 1;
[param, sigma_exact, u_exact, f_exact, uD, gN] = smoothdata2Dmxbd(mu,ld);

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
    [node, elem, misc, misc_boundary] = unitsquaremeshmxbd(N1);
    
    tic;
    %% DOF
    elem2dof_sigma = dofHZ2DP3(elem);
    elem2dof_u = dofDG2D(elem,2);
    DoFtensor = symtensor2HZP3(elem,node);
    
    %% matrices
    C = compliance2(node,elem,elem2dof_sigma,DoFtensor,param);
    B = div2(node,elem,elem2dof_sigma,elem2dof_u,DoFtensor);
    [Ndof_u,Ndof_sigma] = size(B);
    
    %% right hand side and boundary conditions
    F20 = bulkIntegral2(node, elem, elem2dof_u, f_exact);
    [Bn, gn] = elasNeumannboundaryHZ2DP3(node, gN, elem2dof_sigma, misc, misc_boundary);
    fd = elasDirichletboundaryHZ2DP3(node, uD, elem2dof_sigma, misc, misc_boundary);
    N_constraint = size(Bn,1);
    A = [C,B',Bn';...
        B,sparse(Ndof_u,Ndof_u+N_constraint);...
        Bn,sparse(N_constraint,Ndof_u+N_constraint)];
    F = [fd;-F20;gn];
    assembletime(iter) = toc;

    %% solve the equations
    tic;
    SOLU = A\F;
    solvetime(iter) = toc;
    sigma_h = SOLU(1:Ndof_sigma);
    u_h = SOLU(Ndof_sigma+1:Ndof_sigma+Ndof_u);
    
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

