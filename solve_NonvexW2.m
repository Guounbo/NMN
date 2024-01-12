function   [x,w]   =  solve_NonvexW2( y, opts)
global w1;
lambda1 = 1;
 p = 0.1;
J     =   2;
miny   = min(min(abs(y)));
maxy   = max(max(abs(y)));
if maxy == 0
    x     =   zeros( size(y) );
    w    =   ones(size(y));
    return
end
stepy  = (maxy-miny)/1e4;
spany  = (miny:stepy:maxy);
fx     = zeros(1,size(spany,2));
tau    = ones(size(y));
switch lower(opts.string)
    case lower('lp')
        w  = opts.wlp * (y + 1e-16) .^(opts.wlp-1);
        w  = lambda1*w;
        tau   =  (2*(1-p)*w).^(1/(2-p)) + p*w.*(2*(1-p)*w).^((p-1)/(2-p));
    otherwise
        disp('not supported yet');
end
x     =   zeros( size(y) );
i0    =   find( abs(y) > tau );
if length(i0) >= 1
    y0    =   y(i0);
    t     =   abs(y0);
    wt    =   w(i0);
    for  j  =  1 : J
        switch lower(opts.string)
            case lower('lp')
                t    =  abs(y0) - p*wt.*(t).^(p-1);
            otherwise
                disp('not supported yet');
        end
    end
    x(i0)   =  sign(y0).*t;
end
end