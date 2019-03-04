%tvdm limiter function
function f=minmod(a1,a2,a3)
    if (sign(a1)==sign(a2)&&sign(a1)==sign(a3))
        f=min([abs(a1),abs(a2),abs(a3)])*sign(a1);
    else
        f=0;
    end
end