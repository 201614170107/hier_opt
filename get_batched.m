function [Hu,H0,Hd] = get_batched(h,A,B,C,E,D,F) 
% Compute the impulse-response matrices
% Discrete time state space model
% x(k+1) = A x + B u + E d
%      y = C x + D u + F d
% to impulse response matrices
%      y = Hu*u + H0*x0 + Hd*d
[nx,nu] = size(B);
if nargin <5   
    % do not consider process noise    
    M = zeros(nx*h,nx);
    for ix=1:h
        M(1+nx*(ix-1):ix*nx,:)=A^ix;
    end
    
    AB = zeros(h*nx,h*nu);
    for i=1:h
        BCell = repmat({(A^(i-1))*B}, 1, h+1-i);
        ABdiag = blkdiag(BCell{:});
        AB(1+nx*(i-1):end,1:h*nu-(i-1)*nu) = AB(1+nx*(i-1):end,1:h*nu-(i-1)*nu) + ABdiag;
    end
    Hu = C*AB;
    H0 = C*M;
    Hd = zeros(nx,nu);
else 
    % consider process noise
    nd = size(E,2);
    ny = size(C,1); 
    if nargin <= 6 
        D = zeros(ny,nu); 
        F = zeros(ny,nd); 
    end
    
    H0 = zeros((h-1)*ny,nx);
    Hu = zeros((h-1)*ny,h*nu);
    Hd = zeros((h-1)*ny,h*nd);

    % compute first column with all impulse response coefficients
    T = C;
    k1 = 1;
    k2 = ny;
    for k=1:h
        Hu(k1:k2,1:nu) = T*B;
        Hd(k1:k2,1:nd) = T*E;
        T = T*A;
        H0(k1:k2,1:nx) = T;
        k1 = k1+ny;
        k2 = k2+ny;
    end

    % add extra row with direct output contributions for k = 0
    Hu = [D, zeros(size(D,1),size(Hu,2)-size(D,2)); Hu];
    Hd = [F, zeros(size(F,1),size(Hd,2)-size(F,2)); Hd];
    H0 = [C, zeros(size(C,1),size(H0,2)-size(C,2)); H0];
    
    % copy coefficients and fill out remaining columns
    k1row = ny+1;
    k2row = (h+1)*ny;
    k1col = nu+1;
    k2col = nu+nu;
    kk = (h+1)*ny-ny;
    for k=2:h
        Hu(k1row:k2row,k1col:k2col) = Hu(1:kk,1:nu);
        k1row = k1row+ny;
        k1col = k1col+nu;
        k2col = k2col+nu;
        kk = kk-ny;
    end

    k1row = ny+1;
    k2row = (h+1)*ny;
    k1col = nd+1;
    k2col = nd+nd;
    kk = (h+1)*ny-ny;
    for k=2:h
        Hd(k1row:k2row,k1col:k2col) = Hd(1:kk,1:nd);
        k1row = k1row+ny;
        k1col = k1col+nd;
        k2col = k2col+nd;
        kk = kk-ny;
    end
end
