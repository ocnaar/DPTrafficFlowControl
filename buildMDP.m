% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
% Function to build the MDP model (transition and probability matrices) for 
% a dynamic system.
%
% Arguments
% xdgvals: array with number of uniform intervals for gridding each system state.   
% xmin: array of lower bounds for each system state.    
% xmax: array of upper bounds for each system state.
% udgvals: array with number of uniform intervals for gridding each system input.   
% umin: array of lower bounds for each system input.    
% umax: array of upper bounds for each system input.
% Ts: sample time for system discretization. Set Ts=0 for discrete systems.
% f: function with first order dynamics of the system.
% par: structure with parameters of the system.
% wrap: boolean variable to specify if angles in radians need to be wraped. 
%
% Outputs
% SE: matrix of state transitions. 
% SP: matrix with probabilities for state transitions.  

function [SE,SP] = buildMDP(xdgvals,xmin,xmax,udgvals,umin,umax,Ts,f,par,wrap)

% Total number of MDP states
nx = prod(xdgvals);

% Total number of control inputs 
nu = prod(udgvals);

% Number of Neighboring Vertices in a simplex
nv = length(xdgvals)+1;

% uint32 array for storing stochastic state evolution
SE = zeros(nx,nu*nv,'uint32');

% State ID for labeling out of bounds transitions
outbndID = 2^32-1;

% uint8 array for storing likelihood of state transition
SP = zeros(nx,nu*nv,'uint8');

% uint8 value for a probability of one
probone = 2^8-1;

% Arrays for storing next states
if Ts==0
    % Discrete Model
    X = zeros(length(xdgvals),nx);
else
    % Continuous Model
    X = cell(1,nx);
end

% Arrays for storing Khun Triangulation transitions and their probabilities 
x1sidx = zeros(nx,nv);
x1prob = zeros(nx,nv);

% Initialize variable for keeping track of total running time.
tStart = tic;

fprintf('Build MDP Transition Model for %d states and %d inputs.\n',nx,nu);
for j=1:nu
    fprintf('Calculating transitions for input %d...\n',j);
    uj = vidx2vertex(udgvals,sidx2vidx(udgvals,j-1),umin,umax);
    tic
    if Ts==0
        % Discrete Model. 
        parfor i=1:nx
            X(:,i)= f(vidx2vertex(xdgvals,sidx2vidx(xdgvals,i-1),xmin,xmax),uj,par);
        end
        for i=1:nx
            x0 = X(:,i);
            if wrap
                x0wrap = [wrapTo2Pi(x0(1:length(x0)/2)); x0(length(x0)/2+1:end)];
                [outbnd, NN, NNprob] = KhunNN(xdgvals,xmin,xmax,x0wrap);
            else
                [outbnd, NN, NNprob] = KhunNN(xdgvals,xmin,xmax,x0);
            end
            if outbnd
                x1sidx(i,:) = outbndID;
                x1prob(i,:) = probone;
            else
                x1sidx(i,:) = vertex2sidx(xdgvals,xmin,xmax,NN);
                x1prob(i,:) = probone*NNprob;
            end
        end
    else
        % Continuous Model. 
        parfor i=1:nx
            [~,X{i}] = ode45(@(t,x) f(t,x,uj,par),[0 Ts],...
                vidx2vertex(xdgvals,sidx2vidx(xdgvals,i-1),xmin,xmax));
        end
        for i=1:nx
            x0 = X{i}(end,:)';
            if wrap
                x0wrap = [wrapTo2Pi(x0(1:length(x0)/2)); x0(length(x0)/2+1:end)];
                [outbnd, NN, NNprob] = KhunNN(xdgvals,xmin,xmax,x0wrap);
            else
                [outbnd, NN, NNprob] = KhunNN(xdgvals,xmin,xmax,x0);
            end
            if outbnd
                x1sidx(i,:) = outbndID;
                x1prob(i,:) = probone;
            else
                x1sidx(i,:) = vertex2sidx(xdgvals,xmin,xmax,NN);
                x1prob(i,:) = probone*NNprob;
            end
        end
    end
    SE(:,nv*(j-1)+1:nv*j) = x1sidx;
    SP(:,nv*(j-1)+1:nv*j) = x1prob;
    toc
end
tTotal = toc(tStart);
fprintf('\nMDP Transition Model successfully built!\nTotal running time was %d[s] (%d[h]).\n',tTotal,tTotal/3600);