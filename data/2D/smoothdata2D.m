function [param, sigma_exact, u_exact, f_exact] = smoothdata2D(mu, ld)
%SMOOTHDATA2D Summary of this function goes here
%   Detailed explanation goes here

%% problem parameters
param = struct('mu',mu,'ld',ld);

%% exact solution and RHS
syms x y

% exact displacement
u1_exact = exp(x^2)*x*(1-x)^2*y*(1-y)^2; % homogeneous Dirichlet BC
u2_exact = sin(1+x)*x^2*(1-x)*y^2*(1-y);

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


%% transfer to matlab function
u1_exact = matlabFunction(u1_exact);
u2_exact = matlabFunction(u2_exact);

sigma11_exact = matlabFunction(sigma11_exact);
sigma12_exact = matlabFunction(sigma12_exact);
sigma22_exact = matlabFunction(sigma22_exact);

f1_exact = matlabFunction(f1_exact);
f2_exact = matlabFunction(f2_exact);

%% 
sigma_exact = struct('sigma11_exact',sigma11_exact,...
    'sigma12_exact',sigma12_exact,'sigma22_exact',sigma22_exact);
u_exact = struct('u1_exact',u1_exact,...
    'u2_exact',u2_exact);
f_exact = struct('f1_exact',f1_exact,'f2_exact',f2_exact);



end

