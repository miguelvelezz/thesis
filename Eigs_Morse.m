close all
clear all

% This file plots the eigenvalues given by our equation for various
% attraction strengths and verifies the eigenvalues with other numerical
% methods. Also, the pairwise potential energy is global in attraction
% and repulsion. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 51   ;   % Number of agents
L = 100  ;   % Length of interval
D = 1.5  ;   % Range of attraction potential

F_i = linspace(0,1,30); % Attraction strength values that we will iterate through
stability = zeros(1,length(F_i)); % stability of the system at every value of F

% For comparison with EnergyMin code (but maybe not the best parameter values for figure)

% N = 15  ;    % Number of agents
% L = 20  ;    % Length of interval
% D = 10  ;    % Range of attraction potential
% 
% F_i = linspace(0,1,200); % Attraction strength values that we will iterate through
% stability = zeros(1,length(F_i)); % stability of the system at every value of F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting Analytical Eigenvalues %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

titlesize = 16;
axislabelsize = 14;
fontSize = axislabelsize;
ticklabelsize = 10;

%linewidth
lwidth = 1.5;

h =  figure(1);
set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual')
set(h,'PaperPosition',[ 0 0 4 3]);  
set(h,'Position',[ 0 0 4 3]);
get(h,'Position')
hold on
for i = 1:length(F_i)
    eigs = eig_vals(N,L,F_i(i),D);
    stability(i) = sign(min(eigs));
    plot(F_i(i),eigs,'*')
end
yline(0, '-');
set(gca, 'FontSize', ticklabelsize)
xlabel('Attraction strength (F)','FontSize',axislabelsize)
ylabel('Eigenvalues','FontSize',axislabelsize)
title('Stability with Morse potential','FontSize',titlesize)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('eigs_morse');
print(h,[filename, '.eps'],'-depsc') % Look for eps to pdf package in latex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison to Numerical Eigenvalues %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefining parameters (all static)
N = 31     ;
L = 20     ;
F = 0.55   ;
D = 1.5    ;

% Numerical eigenvalues
% num_eig1 uses eig() function / num_eig2 uses fft() function
[num_eig1,num_eig2] = num_eigs(N,L,F,D);

% Analytical eigenvalues with same parameters
analytic_eigs = eig_vals(N,L,F,D);

% The maximum absolute value of the difference between analytical and
% numerical eigenvalues. These values are a good measure of how well my
% analytical calculation agrees with numerical methods
verification1 =  max(abs(analytic_eigs-num_eig1));
verification2 =  max(abs(analytic_eigs-num_eig2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second derivative of Morse potential of attraction and repulsion terms
function secondDerivPotential = Qpp(L,D,F,r) 
    r2 = r-(L*floor(r/L));
    repel = cosh( (L/2) - r2 )/sinh(L/2) ;
    attract = - (F/D)*( cosh( ((L/2) - r2)/D )/sinh(L/(2*D))) ;
    secondDerivPotential = attract + repel;
end

% Constructs Hessian and uses built-in function to compute numerical
% eigenvalues
function [num1,num2] = num_eigs(N,L,F,D)
    delta = L/N ;   % Spacing between particles
    
    % Diagonal value of Hessian
    k = 1:N-1;
    P = Qpp(L,D,F, (k*delta) ); 
    P = sum(P);
    
    % Off diagonal terms
    l     = 0:N-1;
    [r,c] = meshgrid(l,l) ;
    l     = r-c  ;
    S0    =  - Qpp(L,D,F, (l*delta) ) ;
    Z     = ones(N) - diag( ones(N,1) );
    S     = Z.*S0;
    
    H = eye(N)*P + S;   % The hessian 
    h = H(1,:);         % First row vector of H
    
    % Using eig() to compute eigenvalues
    num1  = real(eig(H)); num1(abs(num1)<10^(-10)) = 0; 
    num1  = sort(num1);
    % Using the first row vector h and using fft() to compute eigenvalues
    num2 = real(fft(h)); num2(abs(num2)<10^(-10)) = 0;
    num2 = transpose(num2); num2  = sort(num2);
end

% Operation to calculate eigenvalues using our analytical results  
function analytic_eigenvalues = eig_vals(N,L,F,D)
    analytic_eigenvalues = zeros(N,1);
    M = (N-1)/2 ; 
    delta = L/N ;
    l = 1:M     ;
    for k = 1:N
        eigenvalue = 2 * Qpp(L,D,F, (l*delta)) .* (1 - cos( ((2*pi)/N) * l*k));
        analytic_eigenvalues(k) = sum(eigenvalue); 
    end
    analytic_eigenvalues = sort(analytic_eigenvalues);
end