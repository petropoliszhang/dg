n=[10 20 40 80 160]; %number of elements
xleft=0;
xright=2*pi;
tleft=0;
tright=2*pi;
k=2; %highest degree of polynomials
for i=1:5
    dx=(xright-xleft)/n(i);
    dt=dx/60;
    D=zeros((k+1)*n(i));
    X0=zeros(n(i)*(k+1),1);
    U=zeros(n(i),k+1); %solution matrix

    for j=1:n(i)
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
                first=@(x) basislj(x).*basisgradmj(x); %first part that in(i)volves integral
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
            D((j-1)*(k+1)+1:j*(k+1),(n(i)-1)*(k+1)+1:n(i)*(k+1))=C;
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

    [sol,t]=RK3( D,X0,dt,tleft,tright,n(i)*(k+1));
    evals=2000000;
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
            uh(p)=basis(l,x(p),xj)*sol((j-1)*(k+1)+l+1,t2)+uh(p);   
        end
    end
    
    error(i)=abs((cos(x)*(u-uh)))/evals;
    error2(i)=sum(abs(u-uh))/evals;
end

for i=2:5
    rate(i-1)=log(error(i-1)/error(i))/log(2)
end