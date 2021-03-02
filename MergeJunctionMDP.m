% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Fall 2016
% Perform Value Iteration for the MDP of a Merge Junction in Coogan Samuel 
% and Arcak Murat. "A Benchmark Problem in Transportation Networks". ACM, 
% 2016. 
%
% System States 
% x(1) = occupancy of link 1.
% x(2) = occupancy of link 2.
% x(3) = occupancy of link 3.
%
% System Inputs
% u(1) = uncontrolled number of vehicles arriving on link 1 from upstream.
% u(2) = uncontrolled number of vehicles arriving on link 2.
% u(3) = metered onramp admittance to freeway (control input).

clear;
clc;
close all;

% Merge Junction Parameters
% A discrete time model is used. The length of the links is 1 [mile] and
% the sample time period is 0.5 [min].
par.c = 40;        % [veh/period]
par.v = 0.5;       % [links/period]
par.w = 0.5/3;     % [links/period]
par.xbar = 320;    % [vehicles]
par.beta = 0.75;
par.alpha = 1;
par.alphabar = 5;

% Vector with number of indices per state.
xdgvals = [5 5 5];
% Vector of lower bounds for states.
xmin = zeros(3,1);
% Vector of upper bounds for states.
xmax = par.xbar*ones(3,1);

% Vector with number of indices per input.
udgvals = [1 1 5];
% Vector of lower bounds for inputs.
umin = [40;5;0];
% Vector of upper bounds for inputs.
% umax = 0.25*par.c*ones(3,1);
umax = [40;5;40];

% Number of vertices per simplex in Khun Triangulation
nv = length(xdgvals) + 1;

% Choose a sample time of 0 for discrete models
Ts = 0;

% Build State Evolution Model
[SE,SP] = buildMDP(xdgvals,xmin,xmax,udgvals,umin,umax,Ts,@MergeJunction,par,0);

% Save workspace with transition model 
% save('MergeJuctionMDP01.mat','-v7.3')

%%
% Load previous calculated data
% load('MergeJuctionMDP06.mat')

% Define time horizon based on final simulation time and sample time.  
H = 20;

% Build Stage Cost
fprintf('\nBuilding Stage Cost array based on Total Travel Time...\n')
tic 
[SC,P] = GridTTT(xdgvals,xmin,xmax);
toc

% Find optimal input sequences and value function. 
% First run in laptop took 22.5 minutes. 
% Run in Mac considering xdgvals=[30 30 30] and udgvals=[30 20 20] took 
% 4231.847817[s] (1.1755[h]).
fprintf('\nFind optimal sequence of inputs using value iteration...\n')
tic
[UF,VF] = buildMDPvaluefunctionTTT(SE,SP,SC,P,H,xdgvals,udgvals);
toc

% Choose a sequence of exogenous inputs for Horizon H
D = [umin(1); umin(2)];
D = repmat(D,1,H);
didx = zeros(1,H);
for i=1:H
    didx(i) = state2gsidx(udgvals(1:2),umin(1:2),umax(1:2),D(:,i));
end

% close all
% Create grid of initial states
ng1 = 11;
ng2 = 11;
ng3 = 11;
x1 = linspace(xmin(1),xmax(1),ng1);
x2 = linspace(xmin(2),xmax(2),ng2);
x3 = linspace(xmin(3),xmax(3),ng3);

% Create grid of control inputs
u = linspace(umin(3),umax(3),udgvals(3));

fprintf('\nSimulating system by interpolation of optimal sequence of inputs...\n')

% Figures to draw state evolution and optimal inputs
figure(1)
hold on
figure(2)
hold all
figure(3)
hold all
figure(4)
hold all
figure(5)
hold all
figure(6)
hold all

tic
for l = 1:ng1
    for k = 1:ng2
        for j = 1:ng3
            x0 = [x1(l); x2(k); x3(j)];
            Xt = zeros(H,nv-1);
            Uopt = zeros(H,1);
            % State evolution with Finite State Machine Optimal Control
            % choosing weighted average of optimal inputs for Nearest Vertices
            % in Khun Triangulation.
            for i=1:H
                [outbnd,NV,prob] = KhunNN(xdgvals,xmin,xmax,x0);
                sidxKhun = vertex2sidx(xdgvals,xmin,xmax,NV);
                UFnan = isnan(UF(sidxKhun+1,didx(i)+1,i));
                infea = sum(UFnan);
                if outbnd || infea==nv
                    draw = 0;
                    break
                else
                    usew = find(~UFnan);
                    useKhun = sidxKhun(usew);
                    w = prob(usew)/sum(prob(usew));
                    uinterp = sum(w.*u(UF(useKhun+1,didx(i)+1,i)));
                    exo = vidx2vertex(udgvals(1:2),sidx2vidx(udgvals(1:2),didx(i)),umin(1:2),umax(1:2));
                    x0 = MergeJunction(x0,[exo;uinterp],par);
                    Xt(i,:)=x0';
                    Uopt(i) = uinterp;
                    draw = 1;
                end
            end
            if draw
                t = 0:H;
                figure(1)
                scatter3(x1(l),x2(k),x3(j),'bx')
                figure(2)
                plot(t,[x1(l); Xt(:,1)])
                figure(3)
                plot(t,[x2(k); Xt(:,2)])
                figure(4)
                plot(t,[x3(j); Xt(:,3)])
                figure(5)
                plot(t(1:end-1),Uopt)
                figure(6)
                plot(t,[x1(l)+x2(k)+x3(j); sum(Xt,2)])
                %             V(k,l) = V(k,l)+z;
            else
                figure(1)
                scatter3(x1(l),x2(k),x3(j),'rx')
                %             V(k,l,j) = NaN;
            end
        end
    end
end
toc

figure(1)
grid
xlabel('x_1(0)')
ylabel('x_2(0)')
zlabel('x_3(0)')
title(['Initial States of Links in Merge Junction (states in red are infeasible) d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off
figure(2)
grid
xlabel('k')
ylabel('x_1(k)')
title(['State Evolution for a Merging Junction with a Metered Onramp d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off
figure(3)
grid
xlabel('k')
ylabel('x_2(k)')
title(['State Evolution for a Merging Junction with a Metered Onramp d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off
figure(4)
grid
xlabel('k')
ylabel('x_3(k)')
title(['State Evolution for a Merging Junction with a Metered Onramp d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off
figure(5)
grid
xlabel('k')
ylabel('u(k)')
title(['Optimal Sequence of Controlled Input d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off
figure(6)
grid
xlabel('k')
ylabel('g(k)')
title(['Stage Cost with TTT performance metric d_1=',num2str(umin(1)),' d_2=',num2str(umin(2))])
hold off