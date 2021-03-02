% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
%
% Gives the MDP state index related to the array of natural numbers 
% (including zero) that identify the state values in the grid for each 
% system state.
%
% Arguments
% dgvals: array with number of uniform intervals for gridding each system 
% state.
% vidx: array of natural numbers (including zero) to identify the state 
% values in the grid for each system state.
%
% Output
% sidx: MDP state index.

function sidx = vidx2sidx(dgvals,vidx)

cp = [fliplr(cumprod(fliplr(dgvals))) 1];
sidx = sum(cp(2:end).*vidx);


