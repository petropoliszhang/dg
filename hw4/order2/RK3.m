function [ sol, t ] = RK3( A,X0,dt,t_initial,t_final,order )
    n=(t_final-t_initial)/dt;
    X(1:order,1)=X0;
    t(1)=t_initial;
    ujbar=zeros(order/3);
    ujtilt=zeros(order/3);
    ujtilt2=zeros(order/3);
    for m=2:ceil(n)+1
        k1=(A*X(1:order,m-1))*dt;
        k2=(A*(X(1:order,m-1)+k1/2))*dt;
        k3=(A*(X(1:order,m-1)-k1+2*k2))*dt;
        X(1:order,m)=X(1:order,m-1)+(k1/6+k2*2/3+k3/6);
        %tvdm limiter here "can be improved using matrix product!!!"
        for j=1:order/3
            ujbar(j)=(X(3*j-2,m)+4*X(3*j-1,m)+X(3*j,m))/6;
            ujtilt(j)=X(3*j,m)-ujbar(j);
            ujtilt2(j)=ujbar(j)-X(3*j-2,m);
        end
        X(1,m)=ujbar(1)-minmod(ujtilt2(1),ujbar(2)-ujbar(1),ujbar(2)-ujbar(1));
        X(2,m)=(4*ujbar(1)-minmod(ujtilt(1),ujbar(2)-ujbar(1),ujbar(2)-ujbar(1))+minmod(ujtilt2(1),ujbar(2)-ujbar(1),ujbar(2)-ujbar(1)))/4;
        X(3,m)=ujbar(1)+minmod(ujtilt(1),ujbar(2)-ujbar(1),ujbar(2)-ujbar(1));
        for p=4:3:order-5
            j=floor(p/3)+1;
            X(p,m)=ujbar(j)-minmod(ujtilt2(j),ujbar(j+1)-ujbar(j),ujbar(j)-ujbar(j-1));
            X(p+1,m)=(4*ujbar(j)-minmod(ujtilt(j),ujbar(j+1)-ujbar(j),ujbar(j)-ujbar(j-1))+minmod(ujtilt2(j),ujbar(j+1)-ujbar(j),ujbar(j)-ujbar(j-1)))/4;
            X(p+2,m)=ujbar(j)+minmod(ujtilt(j),ujbar(j+1)-ujbar(j),ujbar(j)-ujbar(j-1));
        end
        k=order/3;
        X(order-2,m)=ujbar(k)-minmod(ujtilt2(k),ujbar(k)-ujbar(k-1),ujbar(k)-ujbar(k-1));
        X(order-1,m)=(4*ujbar(k)-minmod(ujtilt(k),ujbar(k)-ujbar(k-1),ujbar(k)-ujbar(k-1))+minmod(ujtilt2(k),ujbar(k)-ujbar(k-1),ujbar(k)-ujbar(k-1)))/4;
        X(order,m)=ujbar(k)+minmod(ujtilt(k),ujbar(k)-ujbar(k-1),ujbar(k)-ujbar(k-1));
        t(m)=t(m-1)+dt;
    end
    sol=X;
end