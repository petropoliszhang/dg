function f=basisgrad(l,x,xj)
    if l>0
        f=l*(x-xj).^(l-1);
    else
        f=0;
    end
