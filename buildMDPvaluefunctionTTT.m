% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
% Implementation of the one-step-lookahead value iteration method to obtain 
% the value function and optimal inputs for the MDP states. 
%
% Arguments
% SE: matrix of state transitions.
% SP: matrix with probabilities for state transitions.
% SC: array of stage cost.
% P: array of terminal cost.
% H: time horizon.
% xdgvals: array with number of uniform intervals for gridding each system state.
% udgvals: array with number of uniform intervals for gridding each system input.
%
% Outputs
% UF: array of optimal inputs for MDP states. 
% VF: array of value function for MDP states. 

function [UF,VF] = buildMDPvaluefunctionTTT(SE,SP,SC,P,H,xdgvals,udgvals)

% Number of states
n = prod(xdgvals);

% Number of control inputs
m = udgvals(end);

% Number of exogenous inputs
d = prod(udgvals(1:2));

% Number of neighboring vertices in Khun Triangulation
nv = length(xdgvals)+1;

% Array for cost function.
VF = zeros(n,d,H+1);

% Terminal cost.
VF(:,:,end) = repmat(P,1,d);

% uint32 ID for labeling out of bounds transitions
outbndID = 2^32-1;

% uint8 value for a probability of one
probone = 2^8-1;

% Array for optimal inputs.
UF = zeros(n,d,H);

for t = H:-1:1
    for x = 1:n
        for u = 1:d
            % Pick range of SE matrix columns for which the exogenous inputs stay the same.
            a = m*nv*(u-1)+1;
            b = m*nv*u;
            outbnd = sum(SE(x,a:b)==outbndID);
            if outbnd == m*nv
                UF(x,u,t) = NaN;
                VF(x,u,t) = Inf;
            else
                prob = double(reshape(SP(x,a:b),[nv,m]))/probone;
                if outbnd == 0
                    J = reshape(SC(SE(x,a:b)+1)',[nv,m]);
                    ucost = repmat(SC(x),1,m)+sum(prob.*J,1);
                else
                    SEres = reshape(SE(x,a:b),[nv,m]);
                    ucost = zeros(1,m);
                    for k = 1:m
                        if sum(SEres(:,k)==outbndID)>=1
                            ucost(k) = Inf;
                        else
                            ucost(k) = SC(x)+prob(:,k)'*SC(SEres(:,k)+1);
                        end
                    end
                end
                [val,idx] = min(ucost);
                UF(x,u,t) = idx;
                VF(x,u,t) = val;
            end
        end
    end
end