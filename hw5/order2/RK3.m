function [ sol, t,Y ] = RK3( A,X0,dt,t_initial,t_final,order,problem)
    M=1;
    if (problem==1)
        mm=-1;
    else
        mm=0;
    end
    
    n=(t_final-t_initial)/dt;
    X(1:order,1)=X0;
    t(1)=t_initial;
    Y=zeros(order/3,ceil(n)+1);
    for j=1:order/3
            ujbar(j)=(X(3*j-2,1)+4*X(3*j-1,1)+X(3*j,1))/6;
    end
    for p=1:3:order-2
        j=floor((p-1)/3)+1;
        a=2*(X(p+2,1)+X(p,1)-2*X(p+1,1));
        b=4*X(p+1,1)-X(p+2,1)-3*X(p,1);
        if (a==0||(-b)/2/a<=0||(-b)/2/a>=1)
            Mj=max([X(p,1),X(p+2,1)]);
            mj=min([X(p,1),X(p+2,1)]);
        else
            Mj=max([X(p,1),X(p+2,1),-b^2/4/a+X(p,1)]);
            mj=min([X(p,1),X(p+2,1),-b^2/4/a+X(p,1)]);
        end
        theta=min([1,abs((M-ujbar(j))/(Mj-ujbar(j))),abs((mm-ujbar(j))/(mj-ujbar(j)))]);
        if (theta~=1)
            Y(j,1)=1;
        end
        X(p,1)=theta*(X(p,1)-ujbar(j))+ujbar(j);
        X(p+1,1)=theta*(X(p+1,1)-ujbar(j))+ujbar(j);
        X(p+2,1)=theta*(X(p+2,1)-ujbar(j))+ujbar(j);
    end
    for m=2:ceil(n)+1
        k1=(A*X(1:order,m-1))*dt;
        k2=(A*(X(1:order,m-1)+k1/2))*dt;
        k3=(A*(X(1:order,m-1)-k1+2*k2))*dt;
        X(1:order,m)=X(1:order,m-1)+(k1/6+k2*2/3+k3/6);
        %bound-preserving limiter here
        for j=1:order/3
            ujbar(j)=(X(3*j-2,m)+4*X(3*j-1,m)+X(3*j,m))/6;
        end
        for p=1:3:order-2
            j=floor((p-1)/3)+1;
            a=2*(X(p+2,m)+X(p,m)-2*X(p+1,m));
            b=4*X(p+1,m)-X(p+2,m)-3*X(p,m);
            if (a==0||(-b)/2/a<=0||(-b)/2/a>=1)
                Mj=max([X(p,m),X(p+2,m)]);
                mj=min([X(p,m),X(p+2,m)]);
            else
                Mj=max([X(p,m),X(p+2,m),-b^2/4/a+X(p,m)]);
                mj=min([X(p,m),X(p+2,m),-b^2/4/a+X(p,m)]);
            end
            theta=min([1,abs((M-ujbar(j))/(Mj-ujbar(j))),abs((mm-ujbar(j))/(mj-ujbar(j)))]);
            if (theta~=1)
                Y(j,m)=1;
            end
            X(p,m)=theta*(X(p,m)-ujbar(j))+ujbar(j);
            X(p+1,m)=theta*(X(p+1,m)-ujbar(j))+ujbar(j);
            X(p+2,m)=theta*(X(p+2,m)-ujbar(j))+ujbar(j);
        end
        t(m)=t(m-1)+dt;
    end
    sol=X;
end