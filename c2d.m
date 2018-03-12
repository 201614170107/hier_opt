function [Ad,Bd,Ed] = c2d(A,B,E,Ts)
% discretize ss model 
    [nx,nu] = size(B);
    nd = size(E,2);

    dss = expm([A B E; zeros(nu+nd,nx+nu+nd)]*Ts);

    Ad = dss(1:nx,1:nx);
    Bd = dss(1:nx,nx+1:nx+nu);
    Ed = dss(1:nx,nx+nu+1:nx+nu+nd);
end
