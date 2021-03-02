% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
% Function to translate an array of indices vidx to its MDP state.  
%
% Arguments
% dgvals: array with number of uniform intervals for gridding each system state.
% vidx: array of natural numbers (including zero) to identify the state values in the grid for each system state.
% xmin: array with lower bounds for each system state.  
% xmax: array with upper bounds for each system state. 
%
% Outputs
% x: MDP state. 

function x = vidx2vertex(dgvals,vidx,xmin,xmax)

for i=1:length(dgvals)
   grd = linspace(xmin(i),xmax(i),dgvals(i));
   j = vidx(i);
   x(i,1) = grd(j+1);
end

% %% Code test (paste in script outside of function)
% clc
% dgvals = [3 3];
% xmin = [0 -1];
% xmax = [2 1];
% for i=0:prod(dgvals)-1
%    vidx = sidx2vidx(dgvals,i);
%    x = vidx2vertex(dgvals,vidx,xmin,xmax)' 
% end