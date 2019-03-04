function f=basisgrad3(l,x,xj,h)
    xj0=xj-h/2;
    xj1=xj-h/6;
    xj2=xj+h/6;
    xj3=xj+h/2;
    if (l==0)
        f=6./(xj0-xj1)./(xj0-xj2)./(xj0-xj3);
    elseif (l==1)
        f=6./(xj1-xj0)./(xj1-xj2)./(xj1-xj3);
    elseif (l==2)
        f=6./(xj2-xj0)./(xj2-xj1)./(xj2-xj3);
    else
        f=6./(xj3-xj0)./(xj3-xj1)./(xj3-xj2);
    end
end