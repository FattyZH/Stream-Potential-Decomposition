% This function refers to https://github.com/mathworks/Fast-Poisson-Equation-Solver-using-DCT
% DCT solver could solve the Poisson equation with Neumann condition(=0)
function p = poisson_dct(f,dx,dy)
    if(nargin ==1)
        dx = 1;
        dy = 1;
    end

    [nx,ny] = size(f);

    % modified wavenumber
    kx = 0:nx-1;
    ky = 0:ny-1;
    mwx = 2*(cos(pi*kx/nx)-1)/dx^2;
    mwy = 2*(cos(pi*ky/ny)-1)/dy^2;

    fhat = dct(dct(f)')';

    [MWX, MWY] = ndgrid(mwx,mwy);
    phat = fhat./(MWX+MWY);

    % In Neumann condition, there is no solution when phat(1,1)!=0;
    % Here we fix the mean ( with kx=0,ky=0) to be 0
    phat(1,1) = 0;

    % Inverse 2D DCT
    p = idct(idct(phat)')';
end