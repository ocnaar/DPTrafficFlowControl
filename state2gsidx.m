% Octavio Narváez Aroche
% Berkeley Center for Control and Identification 
% Summer 2016
% Finds the index (a natural number, including 0) identifying the nearest 
% MDP state to a system state x.
%
% Arguments 
% dgvals: array with number of uniform intervals for gridding each system state.
% xmin: array of lower bounds for each system state.    
% xmax: array of upper bounds for each system state.
% x: system state. 
%
% Output
% sidx: index of the nearest MDP state to a system state x.

function sidx = state2gsidx(dgvals,xmin,xmax,x)

for i=1:length(x)
    if x(i)<xmin(i) || x(i)>xmax(i)
        outbnd = 1;
        break
    else
        outbnd = 0;
    end
end

if outbnd
    sidx = NaN;
else
    vidx = zeros(1,length(x));
    for i=1:length(x)
        grd = linspace(xmin(i),xmax(i),dgvals(i));
        for j=1:dgvals(i)-1
            if grd(j)<=x(i) && x(i)<=grd(j+1)
                [~,I] = min([x(i)-grd(j),grd(j+1)-x(i)]);
                I = I-2;
                vidx(i) = j+I;
            end
        end
    end
    sidx = vidx2sidx(dgvals,vidx);
end

% Code test (paste in script outside of function)
% dgvals = [3 5];
% xmin = [1; 5];
% xmax = [4; 10];
% for i=0:dgvals(1)-1
%     for j = 0:dgvals(2)-1
%         sidx = state2gsidx(dgvals,xmin,xmax,[1.5+i; 5.5+j])
%     end
% end