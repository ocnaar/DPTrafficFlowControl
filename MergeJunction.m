% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Summer 2016
% Discrete Dynamics for a Merge Junction as described in the Benchmark 
% Problem in Transportation Networks by Coogan Samuel and Arcak Murat. 
% ACM, 2016. 
% States
% x(1): Occupancy of link 1. 
% x(2): Occupancy of link 2.
% x(3): Occupancy of link 3. 
% Exogenous Inputs
% u(1): Uncontrolled number of vehicles arriving on link 1 from upstream. 
% u(2): Uncontrolled number of vehicles arriving on link 2. 
% Control input
% u(3): Metered onramp admittance to freeway.


function xnext = MergeJunction(x,u,par)

% Parameters 
c = par.c;       % [veh/period]
v = par.v;       % [links/period]
w = par.w;       % [links/period]
xbar = par.xbar; % [vehicles]
beta = par.beta;
alpha = par.alpha;
alphabar = par.alphabar;

% Compute demand and supply functions 
d1 = Demand(x(1),c,v);
d2 = Demand(x(2),c,v);
d3 = Demand(x(3),c,v);
s3 = Supply(x(3),w,xbar);

% Vector of next states 
xnext = zeros(3,1);
xnext(1)= x(1)-min([d1,(alpha/beta)*s3])+u(1);
xnext(2)= x(2)-min([d2,alphabar*s3,u(3)])+u(2);
xnext(3)= x(3)-d3+min([beta*d1,alpha*s3])+min([d2,alphabar*s3,u(3)]);