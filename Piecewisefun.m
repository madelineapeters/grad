function f = Piecewisefun(x,t)
    f = piecewise((x>=t)&(x<x.c),P0*normpdf(x-t,m,s),...
    (t>=x)&(x>=t)&(x-t<=x.c),P0*normpdf(x-t,m,s)*exp(-c*(x-x.c)),...
    (t>=x)&(x>=t)&(x-t>x.c),P0*normpdf(x-t,m,s)*exp(-c*t),...
    (x<t)&(x<x.c)&(t-x>=1-x.c),(B^n)*P0*normpdf(x-t,m,s)*exp(-n*c*(1-x.c)),...
    (x<t)&(x<x.c)&(t-x<1-x.c),(B^n)*P0*normpdf(x-t,m,s)*exp(-c*(t-x-(n-1)*x.c)),...
    (x<t)&(x>=x.c)&(t-x>=1-x.c),(B^n)*P0*normpdf(x-t,m,s)*exp(-n*c*(1-x.c))*exp(-c*(x-x.c)),...
    (x<t)&(x>=x.c)&(t-x<1-x.c),(B^n)*P0*normpdf(x-t,m,s)*exp(-c*(t-n*x.c)));
end

