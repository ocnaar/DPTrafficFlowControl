% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Fall 2016
% Function to calculate the stage and terminal costs according to the Total 
% Travel Time (TTT) performance metric for the MDP states of a three link 
% merge junction in transportation networks. 
%
% Arguments
% xdgvals: array with number of uniform intervals for gridding each system state. 
% xmin: vector of lower bounds for system states.
% xmax: vector of upper bounds for system states.
%
% Outputs
% SC: Array of stage cost.
% P: Array of terminal cost. 

function [SC,P] = GridTTT(xdgvals,xmin,xmax)

% Total number of states in the grid
nx = prod(xdgvals);

% Compute Terminal Cost
P = zeros(nx,1);
for i=1:nx
    vidx = sidx2vidx(xdgvals,i-1);
    vertex = vidx2vertex(xdgvals,vidx,xmin,xmax);
    P(i)= sum(vertex);
end

% For the TTT performance metric, the stage and terminal costs are the
% same.
SC = P;