function[Stream,Potential] = SPdec(u,v,dd)
% dd is interval of mesh
if nargin == 2
    dd = 1;
end
% deal with NAN
unan = isnan(u);
vnan = isnan(v);
u(unan) = 0;
v(vnan) = 0;
% calculate Potential function with divergence field and 0 Dirichlet condition
div = divergence(u,v)/dd;
Potential = zeros(size(div)+2);
Potential(2:end-1,2:end-1) = poisson_dst(div,dd,dd);

% Remove Potential flow field
[uk,vk] = gradient(Potential/dd);
u = u-uk(2:end-1,2:end-1);
v = v-vk(2:end-1,2:end-1);

% calculate Boundary conditions of Steam function and curl field
cur = zeros(size(u));
cur(1,:) = cur(1,:) + cur(1,1) + cumtrapz(-v(1,:));
cur(:,1) = cur(:,1) + cur(1,1) + cumtrapz(u(:,1));
cur(end,:) = cur(end,:) + cur(end,1) + cumtrapz(-v(end,:));
cur(:,end) = cur(:,end) + cur(1,end) + cumtrapz(u(:,end));

cur = (-cur - 2*curl(u,v))/dd;

% In Physical, we are used to Counter sign
Stream = -poisson_dst(cur,dd,dd);
Potential = -Potential(2:end-1,2:end-1);

% deal with NAN
Stream(unan&vnan) = nan;
Potential(unan&vnan) = nan;
end
