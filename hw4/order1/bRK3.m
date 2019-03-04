function [ sol, t ,Y] = bRK3( A,X0,dt,t_initial,t_final,order,dx,M)
    n=(t_final-t_initial)/dt;
    X(1:order,1)=X0;
    t(1)=t_initial;
    Y=zeros(order,ceil(n)+1);
    for m=2:ceil(n)+1
        k1=(A*X(1:order,m-1))*dt;
        k2=(A*(X(1:order,m-1)+k1/2))*dt;
        k3=(A*(X(1:order,m-1)-k1+2*k2))*dt;
        X(1:order,m)=X(1:order,m-1)+(k1/6+k2*2/3+k3/6);
        %tvbm limiter here
        if (X(2,m)-X(1,m)~=minmod2(X(2,m)-X(1,m),X(3,m)-X(1,m),X(3,m)-X(1,m),M,dx))
            Y(2,m)=1;
        end
        X(2,m)=X(1,m)+minmod2(X(2,m)-X(1,m),X(3,m)-X(1,m),X(3,m)-X(1,m),M,dx);
        for j=4:2:order-2
            if (X(j,m)-X(j-1,m)~=minmod2(X(j,m)-X(j-1,m),X(j+1,m)-X(j-1,m),X(j-1,m)-X(j-3,m),M,dx))
                Y(j,m)=1;
            end
            X(j,m)=X(j-1,m)+minmod2(X(j,m)-X(j-1,m),X(j+1,m)-X(j-1,m),X(j-1,m)-X(j-3,m),M,dx);
        end
        if (X(order,m)-X(order-1,m)~=minmod2(X(order,m)-X(order-1,m),X(order-1,m)-X(order-3,m),X(order-1,m)-X(order-3,m),M,dx))
            Y(order,m)=1;
        end
        X(order,m)=X(order-1,m)+minmod2(X(order,m)-X(order-1,m),X(order-1,m)-X(order-3,m),X(order-1,m)-X(order-3,m),M,dx);
        t(m)=t(m-1)+dt;
    end
    sol=X;
end