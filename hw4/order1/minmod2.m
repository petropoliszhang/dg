%tvbm minmod function
function f=minmod2(a1,a2,a3,M,dx)
    if (abs(a1)<=M*dx^2)
        f=a1;
    elseif (sign(a1)==sign(a2)&&sign(a1)==sign(a3))
        f=min([abs(a1),abs(a2),abs(a3)])*sign(a1);
    else
        f=0;
    end
end