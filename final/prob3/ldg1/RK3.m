function [ sol, t ] = RK3( A,X0,dt,t_initial,t_final,order )
    n=(t_final-t_initial)/dt;
    X(1:order,1)=X0;
    t(1)=t_initial;
    for m=2:ceil(n)+1
        k1=(A*X(1:order,1))*dt;
        k2=(A*(X(1:order,1)+k1/2))*dt;
        k3=(A*(X(1:order,1)-k1+2*k2))*dt;
        X(1:order,1)=X(1:order,1)+(k1/6+k2*2/3+k3/6);
        t(m)=t(m-1)+dt;
    end
    sol=X;
end