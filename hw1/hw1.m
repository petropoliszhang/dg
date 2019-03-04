n=100; %number of elements
leftend=0;
rightend=1;
h=(rightend-leftend)/n;
k=2; %highest degree of polynomials
uminus=0; %uh_{j-half}
f=@(x) cos(x);
U=zeros(n,k+1); %solution matrix

for j=1:n
    xjmh=leftend+(j-1)*h; %namely, x_{j-half}
    xjph=leftend+j*h; %namely, x_{j+half}
    xj=(j-0.5)*h+leftend;
    A=zeros(k+1); % A in Au_j=b
    b=zeros(k+1,1); % b in Au_j=b
    for m=0:k
        basismj=@(x) basis(m,x,xj);
        basisgradmj=@(x) basisgrad(m,x,xj);
        for l=0:k
            basislj=@(x) basis(l,x,xj);
            first=@(x) basislj(x).*basisgradmj(x); %first part that involves integral
            A(m+1,l+1)=-integral(first,xjmh,xjph)+basislj(xjph).*basismj(xjph);
        end
        second=@(x) f(x).*basismj(x); %second part that involves integral
        b(m+1)=integral(second,xjmh,xjph)+uminus.*basismj(xjmh);
    end
    U(j,:)=A\b;
    uminus=0;
    for l=0:k
        uminus=uminus+U(j,l+1)*basis(l,xjph,xj);
    end
end

q=1000*n; %number of points to use for numerical integral
x=linspace(0,1,q);
uh=zeros(1,q);
u=sin(x);
for p=1:q-1
    j=floor(x(p)*n)+1;
    xj=(j-0.5)*h+leftend;
    for l=0:k
        uh(p)=basis(l,x(p),xj)*U(j,l+1)+uh(p);
    end
end
j=floor(x(p)*n);
xj=(j-0.5)*h+leftend;
for l=0:k
        uh(q)=basis(l,x(q),xj)*U(j,l+1)+uh(q);
end
plot(x,uh)
error1=sum(abs(u-uh))/q %l1 error
error2=sqrt(sum((u-uh).^2)/q) %l2 error
error3=max(u-uh) %linf error