% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
% Function to obtain the array of natural numbers that identifies the values 
% in the grid for each system state from the index related to an MDP state. 
%
% Arguments
% dgvals: array with number of uniform intervals for gridding each system state.
% N: natural index (including zero) for ordered MDP states.  
%
% Output
% vidx: array of natural numbers (including zero) to identify the state 
% values in the grid for each system state.

function vidx = sidx2vidx(dgvals,N)
ndig = numel(dgvals);
cp = fliplr(cumprod(fliplr(dgvals)));
cp = cp(2:end);
vidx = zeros(1,ndig);
for i=1:ndig-1
   ii = floor(N/cp(i));
   vidx(i) = ii;
   N = N - cp(i)*ii;
end
vidx(end) = N;