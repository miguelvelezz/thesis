close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L      = 20      ;  % Length of interval 
F      = 0.25    ;  % Attraction strength
D      = 1.5     ;  % Attraction range
eq_rho = 10      ;  % Average density of swarm (total mass/length of interval)

% Creates a vector of different average densities of swarm (total mass/length of interval)
sigma_max = 200;
sigma_min = 0.01;
sigma_n   = 50;
sigma = linspace(sigma_min,sigma_max,sigma_n);

% Creates a vector of Fourier modes of our perturbation in wave numbers
wave_n = 50;
n = 2*1:wave_n;

[sigma_grid,n_grid] = meshgrid(sigma,n);
k_grid = (2*pi/L)*n_grid;  % Number of radians per distance

B = ((sigma_grid.^2)./(eq_rho^3)).*Rp(L,(sigma_grid./eq_rho)) + ((sigma_grid.^3)./(2*eq_rho^4)).*Rpp(L,(sigma_grid./eq_rho));
C_n = -(8*(D^2)*F*L)./(L^2 + (2*D*pi.*n_grid).^2);

stability = B + (L/4)*C_n;

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
lwidth = 2;

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
title('Stability under Morse repulsion','FontSize',titlesize)
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

filename = sprintf('stabilityBoundary_morse');
print(h,[filename, '.eps'],'-depsc') % Look for eps to pdf package in latex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First derivative of Morse global repulsion function
function diffRep = Rp(L,r) 
    r2 = r-(L*floor(r/L));
    repel = - sinh( (L/2) - r2 )/sinh(L/2) ;
    diffRep = repel;
end

% Second derivative of Morse global repulsion function
function diff2Rep = Rpp(L,r) 
    r2 = r-(L*floor(r/L));
    repel = cosh( (L/2) - r2 )/sinh(L/2) ;
    diff2Rep = repel;
end