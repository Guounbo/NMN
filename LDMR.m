function [Xn] = LDMR(y, A, D, trls, Xinit, Einit, alpha, imgsize)
rho =1.01;
mu = 0.0001;
opts.string='lp';
opts.wlp = 0.1;
eps_abs = 1e-3;
Max_Iter = 100;
p = imgsize(1);
q = imgsize(2);
n = size(A,2);
Xn = Xinit;
En = Einit;
Dn= A*Xn;
Z1n = zeros(prod(imgsize),1);
Z2n = zeros(prod(imgsize),1);
classnum = length(D);
WA = [];
for i=1:length(trls)
    lb = trls(i);
    WA = [WA D(lb)*A(:,i)];
end
M1 = zeros(n, n);
nc = 0;
for i=1:classnum
   pos = find(trls==i);
   M1(pos,pos) = D(i)*A(:,pos)'*A(:,pos);
   nc  = nc + length(pos);
end
M = pinv(alpha/mu*M1+A'*A)*A';

for iter = 1:Max_Iter  
    Z1o = Z1n;
    Z2o = Z2n;
    Xo = Xn;
    % update E
    c = 1 + (2/mu);
    e_tilde = shrinkageW(y - Dn + Z1o/mu,c);
     m1 = reshape( e_tilde, imgsize) ;
       [AU,SU,VU] = svd(m1,'econ');   
        SU = diag(SU);
        [SU] =  solve_NonvexW2( SU, opts);
       En = AU*diag(SU)*VU';
       En = En(:);
    
     % update D

    tempm1=y - En + A*Xo  + (Z1o-Z2o)/mu;
    m2 = reshape(tempm1/2, imgsize); 
        [AU,SU,VU] = svd(m2,'econ');   
        SU = diag(SU);
        [SU] =  solve_NonvexW2( SU, opts);
       Dn = AU*diag(SU)*VU';
       Dn = Dn(:);
    % update X
    g = Dn + Z2o/mu;
    Xn = M*g;  
    % update Z
    Z1n = Z1o + mu*(y- Dn - En);
    Z2n = Z2o + mu*(Dn-A*Xn);
    mu = mu * rho;
    r = A*Xn - y - En;
    s= norm(Dn-A*Xn,'inf')/norm(A*Xn,'inf');
    convergence = (norm(r,2)<eps_abs) && s<eps_abs;

    if convergence
        break;
    end  
       
end

end
function z = shrinkageW(x,c)
z = x./c;
end
