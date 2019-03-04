function f=basis(l,x,xj,h)
xjph=xj+h/2;
if l==0
    f=-(x-xjph)/h*2;
else
    f=(x-xj)/h*2;
end