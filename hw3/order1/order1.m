clear
clc
n=[10 20 40 80 160]; %number of elements
xleft=0;
xright=2*pi;
tleft=0;
tright=2*pi;
k=1; %highest degree of polynomials
A=[2,0;-2,1];
C=[0,-2;0,1];

for i=1:5
    dx=(xright-xleft)/n(i);
    dt=dx/10;
    D=zeros((k+1)*n(i));
    X0=zeros(n(i)*(k+1),1);
    U=zeros(n(i),k+1); %solution matrix
    B=[-4*dx/3, dx/3;dx/3,-dx/3];
    A2=B\A;
    C2=B\C;
    E=zeros(k+1,k+1);
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
            for m=0:k
                E(m+1,l+1)=basis(l,intpoint(m+1),xj,dx);
            end
        end
        S=sin(intpoint);
        X0((j-1)*(k+1)+1:j*(k+1))=E\S;
    end

    [sol,t]=RK3( D,X0,dt,tleft,tright,n(i)*(k+1));
    evals=2000;
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
    
    error(i)=abs((cos(x)*(u-uh)))/evals;
    error2(i)=sum(abs(u-uh))/evals;
end

for i=2:5
    rate(i-1)=log(error2(i-1)/error2(i))/log(2)
end