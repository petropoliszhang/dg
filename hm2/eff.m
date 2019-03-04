n=10; %number of elements
xleft=0;
xright=2*pi;
tleft=0;
tright=2*pi;
dx=(xright-xleft)/n;
dt=dx/4;
k=1; %highest degree of polynomials
D=zeros((k+1)*n);
X0=zeros(n*(k+1),1);
U=zeros(n,k+1); %solution matrix

for j=1:n
    xjmh=xleft+(j-1)*dx; %namely, x_{j-half}
    xjph=xleft+j*dx; %namely, x_{j+half}
    xj=(j-0.5)*dx+xleft;
    A=zeros(k+1); 
    B=zeros(k+1); 
    C=zeros(k+1);
    for m=0:k
        basismj=@(x) basis(m,x,xj);
        basisgradmj=@(x) basisgrad(m,x,xj);
        for l=0:k
            basislj=@(x) basis(l,x,xj);
            basisljm=@(x) basis(l,x,xj-dx);
            first=@(x) basislj(x).*basisgradmj(x); 
            second=@(x) basislj(x).*basismj(x);
            A(m+1,l+1)=-integral(first,xjmh,xjph)+basislj(xjph).*basismj(xjph);
            B(m+1,l+1)=-integral(second,xjmh,xjph);
            C(m+1,l+1)=-basisljm(xjmh)*basismj(xjmh);
        end
    end
    A=B\A;
    C=B\C;
    D((j-1)*(k+1)+1:j*(k+1),(j-1)*(k+1)+1:j*(k+1))=A;
    if (j>1)
        D((j-1)*(k+1)+1:j*(k+1),(j-2)*(k+1)+1:(j-1)*(k+1))=C;
    else
        D((j-1)*(k+1)+1:j*(k+1),(n-1)*(k+1)+1:n*(k+1))=C;
    end
    intpoint=linspace(xjmh,xjph,k+1)';
    E=zeros(k+1,k+1);
    for l=0:k
        basislj=@(x) basis(l,x,xj);
        for m=0:k
            E(m+1,l+1)=basislj(intpoint(m+1));
        end
    end
    S=sin(intpoint);
    X0((j-1)*(k+1)+1:j*(k+1))=E\S;
end

[sol,t]=RK3( D,X0,dt,tleft,tright,n*(k+1));
time=length(t);
x=linspace(0,2*pi,time);
uh=zeros(time-1);
u=zeros(time-1);
cosx=zeros(time-1);

for i=1:time-1
    cosx(:,i)=cos(x(1:time-1)');
end

for t2=1:time-1
    for p=1:time-1
        u(p,t2)=sin(x(p)-t(t2));
        j=floor(x(p)*n/2/pi)+1;
        xj=(j-0.5)*dx+xleft;
        for l=0:k
            uh(p,t2)=basis(l,x(p),xj)*sol((j-1)*(k+1)+l+1,t2)+uh(p,t2);
        end    
    end
end

q=time-1;
error1=sum(sum(abs(u-uh)))/q/q %l1 error
error2=sqrt(sum(sum((u-uh).^2)/q/q)) %l2 error
error3=max(max(u-uh)) %linf error
error4=abs(sum(sum((u-uh).*cosx))/q/q)