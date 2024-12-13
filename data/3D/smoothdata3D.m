function [param, sigma_exact, u_exact, f_exact] = smoothdata3D(mu, ld)
%SMOOTHDATA3D Summary of this function goes here
%   Detailed explanation goes here
%% problem parameters
param = struct('mu',mu,'ld',ld);

%% exact solution and RHS
syms x y z

% exact displacement
u_tmp = x*(1-x)*y*(1-y)*z*(1-z);
u1_exact = 16*u_tmp;
u2_exact = 32*u_tmp;
u3_exact = 64*u_tmp;
clear u_tmp

% exact strain
strain11_exact = diff(u1_exact,x);
strain12_exact = (diff(u1_exact,y)+diff(u2_exact,x))/2;
strain13_exact = (diff(u1_exact,z)+diff(u3_exact,x))/2;
strain22_exact = diff(u2_exact,y);
strain23_exact = (diff(u2_exact,z)+diff(u3_exact,y))/2;
strain33_exact = diff(u3_exact,z);
tracestrain = strain11_exact + strain22_exact + strain33_exact;

% exact stress
sigma11_exact = 2*mu*strain11_exact + ld*tracestrain;
sigma12_exact = 2*mu*strain12_exact;
sigma13_exact = 2*mu*strain13_exact;
sigma22_exact = 2*mu*strain22_exact + ld*tracestrain;
sigma23_exact = 2*mu*strain23_exact;
sigma33_exact = 2*mu*strain33_exact + ld*tracestrain;

% exact load
f1_exact = -diff(sigma11_exact,x) - diff(sigma12_exact,y) - diff(sigma13_exact,z);
f2_exact = -diff(sigma12_exact,x) - diff(sigma22_exact,y) - diff(sigma23_exact,z);
f3_exact = -diff(sigma13_exact,x) - diff(sigma23_exact,y) - diff(sigma33_exact,z);


%% transform symbol to matlab function
u1_exact = matlabFunction(u1_exact);
u2_exact = matlabFunction(u2_exact);
u3_exact = matlabFunction(u3_exact);
sigma11_exact = matlabFunction(sigma11_exact);
sigma12_exact = matlabFunction(sigma12_exact);
sigma13_exact = matlabFunction(sigma13_exact);
sigma22_exact = matlabFunction(sigma22_exact);
sigma23_exact = matlabFunction(sigma23_exact);
sigma33_exact = matlabFunction(sigma33_exact);

sigma_exact = struct('sigma11_exact',sigma11_exact,...
    'sigma12_exact',sigma12_exact,...
    'sigma13_exact',sigma13_exact,...
    'sigma22_exact',sigma22_exact,...
    'sigma23_exact',sigma23_exact,...
    'sigma33_exact',sigma33_exact);
u_exact = struct('u1_exact',u1_exact,...
    'u2_exact',u2_exact,...
    'u3_exact',u3_exact);
f1_exact = matlabFunction(f1_exact);
f2_exact = matlabFunction(f2_exact);
f3_exact = matlabFunction(f3_exact);
f_exact = struct('f1_exact',f1_exact,...
    'f2_exact',f2_exact,...
    'f3_exact',f3_exact);

end

