clear
clc
n=[10 20 40 80]; %number of elements
testnum=length(n);
xleft=0;
xright=2*pi;
tleft=0;
tright=1;
k=1; %highest degree of polynomials
A=zeros(k+1);
C=zeros(k+1);
B=zeros(k+1);
DD=zeros(k+1);
E=zeros(k+1,k+1);

for i=1:testnum
    dx=(xright-xleft)/n(i);
    dt=dx^3/100;
    D=zeros((k+1)*n(i));
    X0=zeros(n(i)*(k+1),1);
    U=zeros(n(i),k+1); %solution matrix
    for j=1:n(i)
        xjmh=xleft+(j-1)*dx; %namely, x_{j-half}
        xjph=xleft+j*dx; %namely, x_{j+half}
        xj=(j-0.5)*dx+xleft;
        for m=0:k
            basismj=@(x) basis(m,x,xj,dx);
            basisgradmj=@(x) basisgrad(m,x,xj,dx);
            for l=0:k
                basislj=@(x) basis(l,x,xj,dx);
                basisgradlj=@(x) basisgrad(l,x,xj,dx);
                basisljm=@(x) basis(l,x,xj-dx,dx);
                basisljp=@(x) basis(l,x,xj+dx,dx);
                first=@(x) basislj(x).*basisgradmj(x); 
                second=@(x) basislj(x).*basismj(x);
                A1(m+1,l+1)=integral(second,xjmh,xjph);
                A2(m+1,l+1)=A1(m+1,l+1);
                A3(m+1,l+1)=A1(m+1,l+1);
                B1(m+1,l+1)=integral(first,xjmh,xjph)+basislj(xjmh).*basismj(xjmh);
                B2(m+1,l+1)=integral(first,xjmh,xjph)-basislj(xjph).*basismj(xjph);
                B3(m+1,l+1)=integral(first,xjmh,xjph)-basislj(xjph).*basismj(xjph);
                C1(m+1,l+1)=-basisljp(xjph)*basismj(xjph);
                C2(m+1,l+1)=basisljm(xjmh)*basismj(xjmh);
                D3(m+1,l+1)=basisljm(xjmh)*basismj(xjmh);
            end
        end
        AA=-(A1\B1)*(A2\C2)*(A3\D3);
        BB=-(A1\B1)*(A2\B2)*(A3\D3)-(A1\B1)*(A2\C2)*(A3\B3)-(A1\C1)*(A2\C2)*(A3\D3)-A3\D3;
        CC=-(A1\B1)*(A2\B2)*(A3\B3)-(A1\C1)*(A2\B2)*(A3\D3)-(A1\C1)*(A2\C2)*(A3\B3)-A3\B3;
        DD=-(A1\C1)*(A2\B2)*(A3\B3);
        D((j-1)*(k+1)+1:j*(k+1),(j-1)*(k+1)+1:j*(k+1))=CC;
        if (j>2&&j<n(i))
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=BB;
            D((j-1)*(k+1)+1:j*(k+1),j*(k+1)+1:(j+1)*(k+1))=DD;
            D((j-1)*(k+1)+1:j*(k+1),(j-3)*(k+1)+1:(j-2)*(k+1))=AA;
        elseif (j==1)
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-1)*(k+1)+1:n(i)*(k+1))=BB;
            D((j-1)*(k+1)+1:j*(k+1),j*(k+1)+1:(j+1)*(k+1))=DD;
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-2)*(k+1)+1:(n(i)-1)*(k+1))=AA;
        elseif (j==n(i))
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=BB;
            D((j-1)*(k+1)+1:j*(k+1),1:(k+1))=DD;
            D((j-1)*(k+1)+1:j*(k+1),(j-3)*(k+1)+1:(j-2)*(k+1))=AA;
        elseif (j==2)
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=BB;
            D((j-1)*(k+1)+1:j*(k+1),j*(k+1)+1:(j+1)*(k+1))=DD;
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-1)*(k+1)+1:n(i)*(k+1))=AA;
        end
        intpoint=linspace(xj,xjph,k+1)';
        S=sin(intpoint);
        X0((j-1)*(k+1)+1:j*(k+1))=S;
    end
    
    [sol,t]=RK3( D,X0,dt,tleft,tright,n(i)*(k+1)); 
    evals=200000;
    x=linspace(0,2*pi,evals);
    uh=zeros(evals,1);
    u=zeros(evals,1);

    t2=length(t);
    for p=1:evals
        u(p)=sin(x(p)-t(t2));
        j=floor(x(p)*n(i)/2/pi)+1;
        if (j==n(i)+1) 
            j=j-1;
        end
        xj=(j-0.5)*dx+xleft;
        for l=0:k
            uh(p)=basis(l,x(p),xj,dx)*sol((j-1)*(k+1)+l+1,t2)+uh(p);   
        end
    end
    
    error2(i)=sum(abs(u-uh))/evals; %moment error
    error3(i)=sqrt(sum((u-uh).^2)/evals); %l2 error
    error4(i)=max(u-uh); %linf error
end

for i=2:testnum
    rate2(i-1)=log(error2(i-1)/error2(i))/log(2); %rate for moment
    rate3(i-1)=log(error3(i-1)/error3(i))/log(2); %rate for l2
    rate4(i-1)=log(error4(i-1)/error4(i))/log(2); %rate for linf
end

for i=3:testnum
    rrate2(i-1)=-log((error2(i-1)-error2(i))/(error2(i-2)-error2(i-1)))/log(2); %rate for moment
    rrate3(i-1)=-log((error3(i-1)-error3(i))/(error3(i-2)-error3(i-1)))/log(2); %rate for l2
    rrate4(i-1)=-log((error4(i-1)-error4(i))/(error4(i-2)-error4(i-1)))/log(2); %rate for linf
end

plot(x',uh,'o')
hold on
plot(x',u)
legend('uh','u')
title('P1 Solution for N=40 at t=1')
hold off