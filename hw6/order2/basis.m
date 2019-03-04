function f=basis(l,x,xj,h)
    xj0=xj-h/2;
    xj1=xj;
    xj2=xj+h/2;
    if (l==0)
        f=(x-xj1).*(x-xj2)./(xj0-xj1)./(xj0-xj2);
    elseif (l==1)
        f=(x-xj0).*(x-xj2)./(xj1-xj0)./(xj1-xj2);
    else
        f=(x-xj0).*(x-xj1)./(xj2-xj0)./(xj2-xj1);
    end
end