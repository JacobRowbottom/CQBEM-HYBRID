function [rho_exact]=ExactRect(t0, mu, omega, c, L, rho_f, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain length L 
%% c=1/|p|, which is speed for Helm and sqrt(omega)/speed for biharmonic
% Quadrature: midpoint rule 
% N = 1000;
% D_omega = omega/10;
% h = D_omega/N; 
% Omega = omega-D_omega/2+h/2:h:omega+D_omega/2-h/2;
% K = Omega/c + 1i*mu/2;
% u_a = h*sum(abs((exp(-1i*K*L)+exp(1i*K*L))./2./K./sin(K*L)).^2.*Omega.^2)*rho_f/c^2/D_omega;
% u_b = (1+exp(-2*mu*L))/(1-exp(-2*mu*L));
%%Rho0 = u_a/u_b;
Rho0=rho_f;
%Rho0=1/(2*c);%*omega


%pause
Xa=X*sec(t0);

% Compute exact solution at X \in [0,L] 
rho_exact = Rho0*(exp(-mu*Xa)+exp(-mu*(2*L-Xa)))./(1-exp(-2*mu*L));
