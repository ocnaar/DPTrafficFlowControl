% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Summer 2016
% Determine the n+1 Nearest Neighbors of state x, relative to the Khun 
% Triangulation of the n-dimensional state space. 
%
% Arguments
% dgvals: array with number of uniform intervals for gridding each system state.
% xmin: array of lower bounds for each system state.    
% xmax: array of upper bounds for each system state.
% x: system state. 
%
% Outputs
% outbnd: boolean value to indicate if the system state is off bounds.  
% NN: array of n+1 Nearest Neighbors of system state x.
% barcoord: array of baricentric coordinates of system state x.

function [outbnd, NN, barcoord] = KhunNN(dgvals,xmin,xmax,x)

% Dimension of state space
n = length(x);

% Array for Nearest Neighbors
NN = zeros(n,n+1);

% Array for n+1 baricentric coordinates
barcoord = zeros(1,n+1);

% Find 2*n coordinates for Neighboring Vertices of x in the grid
coords = zeros(n,2);
for i=1:n
    if x(i)<xmin(i) || xmax(i)<x(i)
        outbnd = 1;
        break
    else
        grd = linspace(xmin(i),xmax(i),dgvals(i));
        for j=1:length(grd)-1
            if grd(j)<=x(i) && x(i)<=grd(j+1)
                coords(i,1) = grd(j);
                coords(i,2) = grd(j+1);
                outbnd = 0;
                break
            else
                outbnd = 0;
            end
        end
    end
end

% Determine Neighboring Vertices if vector x is whitin the grid
if outbnd
    outbnd = 1;
else
    % Binary table for indexing coordinates
    de = 0:2^n-1;
    bi = de2bi(de,n);
    % Array of Neighboring Vertices
    NV = zeros(2^n,n);
    for j=1:n
        for i=1:2^n
            k = bi(i,j)+1;
            NV(i,j) = coords(j,k);
        end
    end
    % Normalize vector x
    xnorm = (x-coords(:,1))./(coords(:,2)-coords(:,1));
    % Sort normalized coordinates and their indeces in descending order
    [sortxnorm,sortidx] = sort(xnorm);
    sortxnorm = flip(sortxnorm);
    sortidx = flip(sortidx)-1;
    % Identify the n+1 vertices of the simplex where point x is located
    NNidx = zeros(1,n+1);
    for i=1:n
       NNidx(i+1) = NNidx(i)+2^sortidx(i);
    end
    NN = NV(NNidx+1,:)';
    % Calculate barycentric coordinates to model likelihood of transition
    % to vertices of the simplex
    barcoord(1) = 1-sortxnorm(1);
    for i=1:n-1
        barcoord(i+1) = sortxnorm(i)-sortxnorm(i+1);
    end
    barcoord(n+1) = sortxnorm(n);
end