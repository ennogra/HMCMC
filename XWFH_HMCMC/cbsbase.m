function B = cbsbase(N,nx)
% N+1: number of the cubic b spline functions with the 1st and last ends at
% non-zero constants. Use N for reflectivity or non-resonant SAXS; N+1 for
% resonant SAXS. Functions are indiced by B_2,B_3,B_4,...,B_(NB-1).
% N+5: total number of knots (including two on each side to complete the
% cubic spline functions)
% N+3: total number of interior knots (the first and last points
% of the range of interest)

% nx: number of x points in the normalized range in [1,NB] <=> [zmin,zmax]
% B: cubic spline bases of dimension length(z) X N+1

NB = N+3;   % total number of interior knots

x = linspace(1,NB,nx);
B = zeros(nx,N+1); 

% n=2 (1st bspline function)
n = 2;
dx = x-n;
ind = (dx>=-2 & dx<-1);
B(ind,n-1) = 6;
ind = (dx>=-1 & dx<0);
B(ind,n-1) = -dx(ind).^3 - 3*dx(ind).^2 - 3*dx(ind)+5;
ind = (dx>=0 & dx<1);
B(ind,n-1) = 2*dx(ind).^3-3*dx(ind).^2-3*dx(ind)+5;
ind = (dx>=1 & dx<2);
B(ind,n-1) = (2-dx(ind)).^3;

for n=3:NB-2
    dx = x-n;
    ind = (dx>=-2 & dx<-1);
    B(ind,n-1) = (2+dx(ind)).^3;
    ind = (dx>=-1 & dx<0);
    B(ind,n-1) = (2+dx(ind)).^3-4*(1+dx(ind)).^3;
    ind = (dx>=0 & dx<1);
    B(ind,n-1) = (2-dx(ind)).^3-4*(1-dx(ind)).^3;
    ind = (dx>=1 & dx<2);
    B(ind,n-1) = (2-dx(ind)).^3;  
end

n = NB-1;
dx = x-n;
ind = (dx>-2 & dx<-1);
B(ind,n-1) = (2+dx(ind)).^3;
ind = (dx>=-1 & dx<0);
B(ind,n-1) = -2*dx(ind).^3-3*dx(ind).^2+3*dx(ind)+5;
ind = (dx>=0 & dx<1);
B(ind,n-1) = dx(ind).^3-3*dx(ind).^2+3*dx(ind)+5;
ind = (dx>=1 & dx<=2);
B(ind,n-1) = 6;


