function f=basisgrad(l,x,xj,h)
    xj0=xj-h/2;
    xj1=xj-h/6;
    xj2=xj+h/6;
    xj3=xj+h/2;
    if (l==0)
        f=((x-xj2).*(x-xj3)+(x-xj1).*(x-xj3)+(x-xj1).*(x-xj2))./(xj0-xj1)./(xj0-xj2)./(xj0-xj3);
    elseif (l==1)
        f=((x-xj2).*(x-xj3)+(x-xj0).*(x-xj3)+(x-xj0).*(x-xj2))./(xj1-xj0)./(xj1-xj2)./(xj1-xj3);
    elseif (l==2)
        f=((x-xj1).*(x-xj3)+(x-xj0).*(x-xj3)+(x-xj0).*(x-xj1))./(xj2-xj0)./(xj2-xj1)./(xj2-xj3);
    else
        f=((x-xj1).*(x-xj2)+(x-xj0).*(x-xj1)+(x-xj0).*(x-xj2))./(xj3-xj0)./(xj3-xj1)./(xj3-xj2);
    end
end