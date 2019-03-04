%first order tvdm/tvbm limiter for discontinuous problem
clear
clc
n=[10 20 40 80 160 320]; %number of elements
testnum=length(n);
xleft=0;
xright=2*pi;
tleft=0;
tright=2*pi;
k=1; %highest degree of polynomials
A=[2,0;-2,1];
C=[0,-2;0,1];

for i=1:testnum
    dx=(xright-xleft)/n(i);
    dt=dx/10;
    D=zeros((k+1)*n(i));
    X0=zeros(n(i)*(k+1),1);
    U=zeros(n(i),k+1); %solution matrix
    B=[-4*dx/3, dx/3;dx/3,-dx/3];
    A2=B\A;
    C2=B\C;
    E=zeros(k+1,k+1);
    S=zeros(k+1,1);
    for j=1:n(i)
        xjmh=xleft+(j-1)*dx; %namely, x_{j-half}
        xjph=xleft+j*dx; %namely, x_{j+half}
        xj=(j-0.5)*dx+xleft;
        D((j-1)*(k+1)+1:j*(k+1),(j-1)*(k+1)+1:j*(k+1))=A2;
        if (j>1)
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=C2;
        else
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-1)*(k+1)+1:n(i)*(k+1))=C2;
        end
        intpoint=linspace(xjmh,xjph,k+1)';
        for l=0:k
            if (intpoint(l+1)<pi/2 || intpoint(l+1)>3*pi/2)
                S(l+1)=0;
            else
                S(l+1)=1;
            end
            for m=0:k
                E(m+1,l+1)=basis(l,intpoint(m+1),xj,dx);
            end
        end
        X0((j-1)*(k+1)+1:j*(k+1))=E\S;
    end
    M=0;
    %[sol,t,Y]=RK3( D,X0,dt,tleft,tright,n(i)*(k+1)); %using tvdm limiter
    [sol,t,Y]=bRK3( D,X0,dt,tleft,tright,n(i)*(k+1),dx,M); %using tvbm limiter
    
    evals=2000;
    x=linspace(0,2*pi,evals);
    uh=zeros(evals,1);
    u=zeros(evals,1);

    t2=length(t);
    for p=1:evals
        ss=x(p)-t(t2)-floor((x(p)-t(t2))/(2*pi))*2*pi;
        if (ss<pi/2 || ss>3*pi/2)
            u(p)=0;
        else
            u(p)=1;
        end
        j=floor(x(p)*n(i)/2/pi)+1;
        if (j==n(i)+1) 
            j=j-1;
        end
        xj=(j-0.5)*dx+xleft;
        for l=0:k
            uh(p)=basis(l,x(p),xj,dx)*sol((j-1)*(k+1)+l+1,t2)+uh(p);   
        end
    end
    
    error(i)=abs((cos(x)*(u-uh)))/evals; %l1 error
    error2(i)=sum(abs(u-uh))/evals; %moment error
    error3(i)=sqrt(sum((u-uh).^2)/evals); %l2 error
    error4(i)=max(u-uh); %linf error
end

for i=2:testnum
    rate(i-1)=log(error(i-1)/error(i))/log(2); %rate for l1;
    rate2(i-1)=log(error2(i-1)/error2(i))/log(2); %rate for moment
    rate3(i-1)=log(error3(i-1)/error3(i))/log(2); %rate for l2
    rate4(i-1)=log(error4(i-1)/error4(i))/log(2); %rate for linf
end

%below are just for plotting graph
uh2=zeros(evals,length(t));
for tt=1:length(t)
    for p=1:evals
        j=floor(x(p)*n(i)/2/pi)+1;
        if (j==n(i)+1) 
            j=j-1;
        end
        xj=(j-0.5)*dx+xleft;
        for l=0:k
            uh2(p,tt)=basis(l,x(p),xj,dx)*sol((j-1)*(k+1)+l+1,tt)+uh2(p,tt);   
        end
    end
end

x2=linspace(0,2*pi,n(i));
Y=Y(2:2:n(i)*(k+1),:);
q=1;
uh3=zeros(n(i),length(t));
for tt=1:length(t)
    for p=1:n(i)
        j=floor(x2(p)*n(i)/2/pi)+1;
        if (j==n(i)+1) 
            j=j-1;
        end
        xj=(j-0.5)*dx+xleft;
        for l=0:k
            uh3(p,tt)=basis(l,x2(p),xj,dx)*sol((j-1)*(k+1)+l+1,tt)+uh3(p,tt);   
        end
        if(Y(p,tt)==1)
            x3(q)=x2(p);
            y3(q)=t(tt);
            z3(q)=uh3(p,tt);
            q=q+1;
        end
    end
end

[X T]=meshgrid(x,t);
mesh(X',T',uh2)
title('Solution using TVD Limiter')
hold on
plot3(x3,y3,z3,'r.')
hold off