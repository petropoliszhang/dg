function [ sol, t,Y] = RK3( A,X0,dt,t_initial,t_final,order,problem)
    M=1;
    if (problem==1)
        mm=-1;
    else
        mm=0;
    end
    n=(t_final-t_initial)/dt;
    X(1:order,1)=X0;
    t(1)=t_initial;
    Y=zeros(order/2,ceil(n)+1);
    for p=1:2:order-1
        j=floor((p-1)/2)+1;
        Mj=max(2*X(p,1)-X(p+1,1),X(p+1,1));
        mj=min(2*X(p,1)-X(p+1,1),X(p+1,1));
        theta=min([1,abs((M-X(p,1))/(Mj-X(p,1))),abs((mm-X(p,1))/(mj-X(p,1)))]);
        if (theta~=1)
            Y(j,1)=1;
        end
        X(p,1)=theta*(X(p,1)-X(p,1))+X(p,1);
        X(p+1,1)=theta*(X(p+1,1)-X(p,1))+X(p,1);
    end
    for m=2:ceil(n)+1
        k1=(A*X(1:order,m-1))*dt;
        k2=(A*(X(1:order,m-1)+k1/2))*dt;
        k3=(A*(X(1:order,m-1)-k1+2*k2))*dt;
        X(1:order,m)=X(1:order,m-1)+(k1/6+k2*2/3+k3/6);
        %bound-preserving limiter here
        for p=1:2:order-1
            j=floor((p-1)/2)+1;
            Mj=max(2*X(p,m)-X(p+1,m),X(p+1,m));
            mj=min(2*X(p,m)-X(p+1,m),X(p+1,m));
            theta=min([1,abs((M-X(p,m))/(Mj-X(p,m))),abs((mm-X(p,m))/(mj-X(p,m)))]);
            if (theta~=1)
                Y(j,m)=1;
            end
            X(p,m)=theta*(X(p,m)-X(p,m))+X(p,m);
            X(p+1,m)=theta*(X(p+1,m)-X(p,m))+X(p,m);
        end
        t(m)=t(m-1)+dt;
    end
    sol=X;
end