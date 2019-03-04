N = 80;
kwave = 4;
sigma = 0.49;
Tend = 1;
nu = 0.1;

%
h = 1/N;
k = sigma*h^2/abs(nu);

%
x = [h:h:1-h]';
%u = sin(kwave*pi*x);
 u = rand(N-1,1);
time = 0;

%
A = diffmat(N-1);

%
while (time<Tend)
   u = u + nu*k/h^2*A*u;
   time = time+k

   plot(x,u);   axis([0 1 -1 1]);
   drawnow;
   
end;

uexact = sin(kwave*pi*x)*exp(-nu*(pi*kwave)^2*time);

error = sqrt(h*sum( (u-uexact).^2 ))

