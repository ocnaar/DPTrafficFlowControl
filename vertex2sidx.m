% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
%
% Provides the MDP state indices (natural numbers, including 0) identifying 
% the vertices of a simplex in the Khun Triangulation of the state space.
%
% Arguments 
% dgvals: array with number of uniform intervals for gridding each system state.
% xmin: array of lower bounds for each system state.    
% xmax: array of upper bounds for each system state.
% NN: array with vertices of a simplex. 
%
% Output
% sidx: array of MDP state indices for each of the vertices in NN.

function sidx = vertex2sidx(dgvals,xmin,xmax,NN)

% Number of states and vertices
[ns,nv] = size(NN);

% Array for storing sidx values
sidx = zeros(1,nv);

% Array for storing vidx values
vidx = zeros(ns,nv);

for i=1:ns
    grd = linspace(xmin(i),xmax(i),dgvals(i));
    for j=1:nv
        for k=1:dgvals(i)
            if NN(i,j)==grd(k)
                vidx(i,j) = k-1;
                break
            end
        end
    end
end

for i=1:nv
    sidx(i) = vidx2sidx(dgvals,vidx(:,i)');
end

% %% Code test (paste in script outside of function)
% clc
% dgvals = [40 40 40 40];
% xmin = [0 -1 10 -5];
% xmax = [2 1 20 10];
% x = zeros(length(xmin),prod(dgvals));
% 
% for i=0:prod(dgvals)-1
%    vidx = sidx2vidx(dgvals,i);
%    x(:,i+1) = vidx2vertex(dgvals,vidx,xmin,xmax);
% end
% 
% tic
% for i=1:1
%     sidx = vertex2sidx(dgvals,xmin,xmax,x);
% end
% toc