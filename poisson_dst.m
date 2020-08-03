% DST solver could solve the Poisson equation with Dirichlet condition(default 0)
function p = poisson_dst(f,dx,dy)
    if(nargin ==1)
        dx = 1;
        dy = 1;
    end

    [nx,ny] = size(f);

    % modified wavenumber
    kx = 1:nx;
    ky = 1:ny;
    mwx = 2*(cos(pi*kx/(nx+1))-1)/dx^2;
    mwy = 2*(cos(pi*ky/(ny+1))-1)/dy^2;

    % 2D DST
    fhat = dst(dst(f)')';

    [MWX, MWY] = ndgrid(mwx,mwy);
    phat = fhat./(MWX+MWY);

    % Inverse 2D DST
    p = idst(idst(phat)')';
end