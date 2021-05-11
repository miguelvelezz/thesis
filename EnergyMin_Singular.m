close all
clear all
rng(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hysteresis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Varying F in the forward direction
%F_f = [0.01:0.01:0.24 0.25:0.005:0.27 0.28:0.01:1]; F_f = F_f(:); % Transitions between 0.2601-0.2602
F_f = [1:1:300]; F_f = F_f(:);

% Varying F in the reverse direction
F_r = [0.01:0.01:0.16 0.19:0.005:0.24 0.25:0.01:1];                % Transitions between 0.221-0.220 but with noise
F_r = F_r(end:-1:1); F_r = F_r(:);

% System analysis with F forward
[supp_f, max_grad_f, eigenvalues_f] = system_analysis(F_f);
e_f = min(eigenvalues_f); min_eig_f = e_f(:); max_grad_f(max_grad_f>0.005) = NaN; supp_f(supp_f>1000) = NaN;

% System analysis with F reverse
%[supp_r, max_grad_r, eigenvalues_r] = system_analysis(F_r);
%e_r = min(eigenvalues_r); min_eig_r = e_r(:); max_grad_r(max_grad_r>0.005) = NaN; supp_r(supp_r>1000) = NaN;

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

h = figure(2);
hold on
plot(F_f,supp_f,'bo')
%plot(F_r,supp_r,'r+')
xlabel('Attraction strength (F)','FontSize',axislabelsize) 
ylabel('Support','FontSize',axislabelsize)
title('Swarm size vs attraction strength','FontSize',titlesize)
hold off

figure(3)
hold on
plot(F_f,max_grad_f,'bo')
%plot(F_r,max_grad_r,'r+')
xlabel('F','FontSize',18) 
ylabel('Max gradient value','FontSize',18)
hold off

% figure(4)
% hold on
% plot(F_f,min_eig_f,'bo')
% plot(F_r,min_eig_r,'r+')
% xlabel('F','FontSize',18) 
% ylabel('Min eigenvalue','FontSize',18)
% hold off

set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual')
set(h,'PaperPosition',[ 0 0 7 5]); % (0,0) - bottom left location; (4,3) 4 wide, 3 tall 
set(h,'Position',[ 0 0 7 5]);
get(h,'Position')

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename = sprintf('hysteresis_singular');
% print(h,[filename, '.eps'],'-depsc') % Look for eps to pdf package in latex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective function - Outputs energy,gradient,hessian of a particle
% configuration x

function [Energy,gradE,H] =W(n,F,D,L,L_r,x)
    epsilon = 0.0001;
    [Xj, Xi] = meshgrid(x, x);
    Xij = Xj - Xi;
%   Energy
    e = R(L,L_r,Xij) + A(L,D,F,Xij);  % Repulsion + Attraction
    e(1:1+size(e,1):end) = 0;         % diagonal = 0; no self interaction
    eout = e;                         % Matrix of pairwise interaction energies
    Energy = 1/2 * sum(sum(eout)) + Phi(x,L,epsilon); % superposition of interaction energies + confining potential

%   Gradient 
    epout = Rp(L,L_r,Xij) + Ap(L,D,F,Xij);
    epout(1:1+size(epout,1):end) = 0;      % diagonal = 0
    dEdx= sum(epout,1);                    % Gradient of social forces
    gradE = dEdx(:) + Phip(x,L,epsilon);   % Gradient of social forces + gradient of confining potential
   
%   Hessian 
    eppout = Rpp(L,L_r,Xij) + App(L,D,F,Xij);
    eppout(1:1+size(eppout,1):end) = 0;    % diagonal = 0
    Z = ones(n) - diag( ones(n,1) );       % Matrix of ones except for zeros on diagonal
    S = Z.*eppout;                         % S is a matrix of off diagonal terms
    P = diag(sum(S));                      % Diagonal term
    H = P-S + diag(2*epsilon*ones(n,1));   % Both matrices together with the contribution of confining potential
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confining potential (Mimics center of mass constraint)

function Confining = Phi(x,L,epsilon)
    x = x - L/2;
    Confining = epsilon * sum(x.^2);
end

function diffConfining = Phip(x,L,epsilon)
    x = x - L/2;
    diffConfining = 2 * epsilon * x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Singular repulsion

function Rep = R(L,L_r,r)
    r2    = r-(L*floor(r/L));
    r3    = abs(r2+(r2>L/2).*(-L));
    p     = 2;
    repel = ( (L_r./r3).^p + (p/(1-p)).*(L_r./r3).^(p-1) - 1/(p-1) ) .* (1/L_r)^p .* heaviside(L_r-r3) .*  heaviside(L_r+r3);
    repel(abs(r3)<eps) = 0;
    Rep   = repel;
end

function diffRep = Rp(L,L_r,r)
    r2    = r-(L*floor(r/L));
    r3    = r2+(r2>L/2).*(-L);   
    p     = 2;
    repel = sign(r3).* (-p).* (L_r./abs(r3)).^(p+1) * (1/L_r)^(p+1) .* (1- (abs(r3)./L_r)) .* heaviside(L_r-r3) .*  heaviside(L_r+r3);
    repel(abs(r3)<eps) = 0;
    diffRep = repel;
end

function diff2Rep = Rpp(L,L_r,r)
    r2    = r-(L*floor(r/L));
    r3    = abs(r2+(r2>L/2).*(-L));
    p     = 2;
    repel = (1/L_r)^(p+2) .* ( (p+1)*p*(L_r./r3).^(p+2) .* (1- (r3./L_r)) + p*(L_r./r3).^(p+1) ) .* heaviside(L_r-r3) .*  heaviside(L_r+r3);
    repel(abs(r3)<eps) = 0;
    diff2Rep = repel;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Morse Attraction

function Att = A(L,D,F,r)
    r2 = r-(L*floor(r/L));
    attract = - F*D*( cosh( ((L/2) - r2)/D )/sinh(L/(2*D))) ;
    attract(abs(r)<eps) = 0;
    Att = attract;
end

function diffAtt = Ap(L,D,F,r)
    r2 = r-(L*floor(r/L));
    attract = F*( sinh( ((L/2) - r2)/D )/sinh(L/(2*D))) ;
    attract(abs(r)<eps) = 0;
    diffAtt = attract;
end

function diff2Att = App(L,D,F,r) 
    r2 = r-(L*floor(r/L));
    attract = - (F/D)*( cosh( ((L/2) - r2)/D )/sinh(L/(2*D))) ;
    attract(abs(r)<eps) = 0;
    diff2Att = attract;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System analysis varies a parameter and computes information about the swarm equilibrium  
% supp - the support of our swarm density, max_grad - the maximum value of
% the gradient lets us know if the system has converged to an equilibrium,
% eigenvalues provides stability information of the swarm.  

function [supp, max_grad, eigenvalues] = system_analysis(lambda)
    n = 15     ;  
    D = 10     ; 
    L = 20     ;
    L_r = 6    ;
    supp        = zeros(length(lambda),1);
    max_grad    = zeros(length(lambda),1);
    eigenvalues = zeros(n,length(lambda));

    x0 = linspace(0,L,n+1);
    x0 = x0(1:n);
    x0 = x0-mean(x0)+(L/2);
    x  = x0' + (rand(size(x0'))-0.5)*0.001;
    
    lb = zeros(size(x0(:)));
    ub = L*ones(size(x0(:)));
    
    options = optimoptions('fmincon', 'Algorithm','trust-region-reflective',...
        'SpecifyObjectiveGradient',true,'HessianFcn','objective', 'MaxIterations', 5e4 ,...
        'TolX', 1e-16, 'TolFun', 1e-14, 'MaxFunctionEvaluations', 5e4);
    
    count_red = 0;
    figure(1)
    hold on
    for i = 1:length(lambda)
        xs = fmincon(@(input) W(n,lambda(i),D,L,L_r,input),x,[],[],[],[],lb,ub,[],options);  % Optimization routine, returns minimum energy configuration
        xs = sort(xs);
        if i == 1
            xold = xs;
        end
        [Energy,gradE,H] = W(n,lambda(i),D,L,L_r,xs);   % Retreiving energy, gradient, and hessian information from equilibrium xs
        Energy
        max_grad(i) = max(abs(gradE));                  % Save the maximum magnitude gradient value to check if near zero
        eigenvalues(:,i) = eig(H);                      % Save eigenvalues to evaluate stability
        supp(i) = abs(xs(1)-xs(end)) + (max(abs(gradE))>0.005)*1000;        % Save the support of the equilibrium configuration
        xaxis = lambda(i)*ones(size(xs));               % x-axis of plot of equilibrium configuration vs parameter values (F)
        if max(abs(gradE))<0.005
            %plot(xaxis(:),xs,'b+') 
        else
            count_red = count_red + 1;
            %plot(xaxis(:),xs,'r+')
        end
        plot(xaxis(:),xs,'b+')
        yline(L/2, '-');                                 % Plots horizontal line at center of mass of swarm
        xlabel('Attraction strength (F)','FontSize',18) 
        ylabel('Particle positions','FontSize',18)
        title('Equilibrium positions','FontSize',16)
        x = 2*xs - xold + (rand(size(x0'))-0.5)*0.00;  % Secant method
        xold = xs;
    end
    count_red
end