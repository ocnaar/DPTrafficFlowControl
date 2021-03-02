% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Summer 2016
% Supply function for a triangular fundamental diagram in transportation 
% networks. From Coogan Samuel and Arcak Murat. A Benchmark Problem in 
% Transportation Networks. ACM, 2016. 
%
% Arguments
% x: occupancy of the link [vehicles].
% w: congestion-wave speed [links/period]. 
% xbar: jam occupancy of the link [vehicles]. 
%
% Output
% s: supply of a link [vehicles/period]. 

function s = Supply(x,w,xbar)
s = max([w*(xbar-x),0]);