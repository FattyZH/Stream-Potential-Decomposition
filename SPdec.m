function[Stream,Potential] = SPdec(u,v)

unan = isnan(u);
vnan = isnan(v);
u(unan) = 0;
v(vnan) = 0;

div = divergence(u,v);
Potential = zeros(size(div)+2);
Potential(2:end-1,2:end-1) = poisson_dst(div);
[uk,vk] = gradient(Potential);
u = u-uk(2:end-1,2:end-1);
v = v-vk(2:end-1,2:end-1);


ul = -(u(1:end-1,1)+u(2:end,1))/2;ur = (u(1:end-1,end)+u(2:end,end))/2;
vb = -(v(1,1:end-1)+v(1,2:end))/2;vu = (v(end,1:end-1)+v(end,2:end))/2;
vm = mean([ul(:);ur(:);vb(:);vu(:)]);
ul = ul-vm; ur = ur-vm; vb = vb-vm;vu = vu-vm;
cur = zeros(size(u));
for i = 1:size(u,1)-1
    cur(i+1,1) = cur(i,1)-ul(i);
end
for i = 1:size(u,2)-1
    cur(1,i+1) = cur(1,i)+vb(i);
end
for i = 1:size(u,1)-1
    cur(i+1,end) = cur(i,end)+ur(i);
end
for i = 1:size(u,2)-2
    cur(end,i+1) = cur(end,i)-vu(i);
end
cur(1,1) = 2*cur(1,1);cur(1,end) = 2*cur(1,end);cur(end,1) = 2*cur(end,1);cur(end,end) = 2*cur(end,end);
cur = -cur - 2*curl(u,v);

Stream = poisson_dst(cur);
Potential = Potential(2:end-1,2:end-1);

Stream(unan&vnan) = nan;
Potential(unan&vnan) = nan;
end
