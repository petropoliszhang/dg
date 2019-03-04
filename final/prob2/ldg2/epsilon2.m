clear
clc
n=[10 20 40 80 160 320]; %number of elements
testnum=length(n);
xleft=0;
xright=2*pi;
tleft=0;
tright=1;
k=2; %highest degree of polynomials
A=zeros(k+1);
C=zeros(k+1);
B=zeros(k+1);
DD=zeros(k+1);
E=zeros(k+1,k+1);
epsilon=0.01;

for i=1:testnum
    dx=(xright-xleft)/n(i);
    dt=dx^2/70;
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
                A(m+1,l+1)=integral(first,xjmh,xjph)-basislj(xjph).*basismj(xjph);
                B(m+1,l+1)=-integral(second,xjmh,xjph);
                C(m+1,l+1)=basisljm(xjmh)*basismj(xjmh);
                A4(m+1,l+1)=epsilon*integral(first,xjmh,xjph)+integral(second,xjmh,xjph)+epsilon*basislj(xjmh).*basismj(xjmh);
                D4(m+1,l+1)=-basisljp(xjph)*basismj(xjph)*epsilon;
            end
        end
        A2=B\A;
        C2=B\C;
        DD2=B\DD;
        A4=B\A4;
        D4=B\D4;
        B3=A4*C2;
        C3=A4*A2+D4*C2;
        D3=D4*A2;
        D((j-1)*(k+1)+1:j*(k+1),(j-1)*(k+1)+1:j*(k+1))=C3;
        if (j>1&&j<n(i))
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=B3;
            D((j-1)*(k+1)+1:j*(k+1),j*(k+1)+1:(j+1)*(k+1))=D3;
        elseif (j==1)
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-1)*(k+1)+1:n(i)*(k+1))=B3;
            D((j-1)*(k+1)+1:j*(k+1),j*(k+1)+1:(j+1)*(k+1))=D3;
        elseif (j==n(i))
            D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=B3;
            D((j-1)*(k+1)+1:j*(k+1),1:(k+1))=D3;
        end
        intpoint=linspace(xjmh,xjph,k+1)';
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
        u(p)=sin(x(p)-t(t2))*exp(-t(t2)*epsilon);
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
hold off