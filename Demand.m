% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Summer 2016
% Demand function for a triangular fundamental diagram in transportation 
% networks. From Coogan Samuel and Arcak Murat. A Benchmark Problem in 
% Transportation Networks. ACM, 2016. 
%
% Arguments
% x: occupancy of the link [vehicles].
% c: capacity of the link [vehicles/period]. 
% v: free-flow speed [links/period]. 
%
% Output
% d: demand of the link [vehicles/period]. 

function d = Demand(x,c,v)
d = min([c,v*x]);