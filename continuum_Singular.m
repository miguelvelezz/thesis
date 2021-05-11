close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L      = 100     ;  % Length of interval
L_r    = 20      ;  % Repulsion range
F      = 0.25    ;  % Attraction strength
D      = 3.00    ;  % Attraction range 
eq_rho = 10      ;  % Equilibrium constant density value (rho_bar in thesis)

% Creates a vector of different average densities of swarm (total mass/length of interval)
sigma_max = 10;
sigma_min = 0.01;
sigma_n   = 70;
sigma = linspace(sigma_min,sigma_max,sigma_n);

% Creates a vector of Fourier modes of our perturbation in wave numbers
wave_n = 50;
n = 2*1:wave_n;

[sigma_grid,n_grid] = meshgrid(sigma,n);
k_grid = (2*pi/L)*n_grid;  % Number of radians per distance

B = ((sigma_grid.^2)./(eq_rho^3)).*Rp(L,L_r,(sigma_grid./eq_rho)) + ((sigma_grid.^3)./(2*eq_rho^4)).*Rpp(L,L_r,(sigma_grid./eq_rho));
C_n = -(8*(D^2)*F*L)./(L^2 + (2*D*pi.*n_grid).^2);

stability = B + (L/4)*C_n;  % k^2 proportional to rho 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
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
set(h,'PaperPosition',[ 0 0 7 5]);  
set(h,'Position',[ 0 0 7 5]);
get(h,'Position')
hold on
surf(sigma_grid,k_grid,sign(stability))
%colorbar
set(gca, 'FontSize', ticklabelsize)
xlabel('Average density (M/L)','FontSize',axislabelsize)
ylabel('Wavenumber (k)','FontSize',axislabelsize)
title('Stability with singular repulsion','FontSize',titlesize)
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

filename = sprintf('stabilityBoundary_singular');
print(h,[filename, '.eps'],'-depsc') % Look for eps to pdf package in latex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First derivative of singular local repulsion function
function diffRep = Rp(L,L_r,r)
    r2    = r-(L*floor(r/L));
    r3    = r2+(r2>L/2).*(-L);   
    p     = 2;
    repel = sign(r3).* (-p).* (L_r./abs(r3)).^(p+1) * (1/L_r)^(p+1) .* (1- (abs(r3)./L_r)) .* heaviside(L_r-r3) .*  heaviside(L_r+r3);
    repel(abs(r3)<eps) = 0;
    diffRep = repel;
end

% Second derivative of singular local repulsion function
function diff2Rep = Rpp(L,L_r,r)
    r2    = r-(L*floor(r/L));
    r3    = abs(r2+(r2>L/2).*(-L));
    p     = 2;
    repel = (1/L_r)^(p+2) .* ( (p+1)*p*(L_r./r3).^(p+2) .* (1- (r3./L_r)) + p*(L_r./r3).^(p+1) ) .* heaviside(L_r-r3) .*  heaviside(L_r+r3);
    repel(abs(r3)<eps) = 0;
    diff2Rep = repel;
end