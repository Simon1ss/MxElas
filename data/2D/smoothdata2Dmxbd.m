function [param, sigma_exact, u_exact, f_exact, uD ,gN] = smoothdata2Dmxbd(mu, ld)
%SMOOTHDATA2DMXBD Summary of this function goes here
%   Detailed explanation goes here
%% artificial problem
% domain \Omega = [0,1]^2
% Neumann boundary: x=0 | x=1
% Dirichlet boundary: y=0 | y=1

%% problem parameters
param = struct('mu',mu,'ld',ld);

%% exact solution
syms x y
u1_exact = sin(1+x^2)*exp(y);
u2_exact = exp(x)*cos(y^2);


% (sigma-u): P3-P4 !!!!!!!!!! Focused on this problem
% u1 = x^2*y^2;
% u2 = x^3*y;

% (sigma-u): P2-P3
%u1 = randn*x^2*y+x*y^2;
%u2 = x^3+y^3;
% u1 = x^2*y+x*y^2;
% u2 = x^3+y^3;



% u1 = y*(1-y)*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);
% u2 = y*(1-y)*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);

% u1 = y*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);
% u2 = y*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);
% 


% u1 = (1-y)*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);
% u2 = (1-y)*(randn*x^3+randn*x^2*y+randn*y^3+randn*x*y^2+...
%     randn*x^2+randn*x*y+randn*y^2+randn*x+randn*y+randn);
% 



%% solid strain
strain11_exact = diff(u1_exact,x);
strain12_exact = (diff(u1_exact,y)+diff(u2_exact,x))/2;
strain22_exact = diff(u2_exact,y);
tracestrain = strain11_exact+strain22_exact;

%% solid stress
sigma11_exact = 2*mu*strain11_exact + ld*tracestrain;
sigma12_exact = 2*mu*strain12_exact;
sigma22_exact = 2*mu*strain22_exact + ld*tracestrain;

%% solid load
f1_exact = - diff(sigma11_exact,x)-diff(sigma12_exact,y);
f2_exact = - diff(sigma12_exact,x)-diff(sigma22_exact,y);


%% boundary conditions
% Dirichlet boundary
% displacement
u1_exact = matlabFunction(u1_exact); % function of (x,y)
u2_exact = matlabFunction(u2_exact);

% Neumann boundary
nx = @(x) 2*x-1;
sigma11_exact = matlabFunction(sigma11_exact);
sigma12_exact = matlabFunction(sigma12_exact);
sigma22_exact = matlabFunction(sigma22_exact);
% Neumann boundary data
gN1_exact = @(x,y) nx(x).*sigma11_exact(x,y);
gN2_exact = @(x,y) nx(x).*sigma12_exact(x,y);

%% transfer to matlab function
f1_exact = matlabFunction(f1_exact);
f2_exact = matlabFunction(f2_exact);

%% restore in structures
sigma_exact = struct('sigma11_exact',sigma11_exact,...
    'sigma12_exact',sigma12_exact,'sigma22_exact',sigma22_exact);
u_exact = struct('u1_exact',u1_exact,'u2_exact',u2_exact);
f_exact = struct('f1_exact',f1_exact,'f2_exact',f2_exact);
uD = struct('uD1_exact',u1_exact,'uD2_exact',u2_exact);
gN = struct('gN1_exact',gN1_exact,'gN2_exact',gN2_exact);






end

