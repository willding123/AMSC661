% Solve PDE from stationary heat distribution 

fsz = 20;

% generate meshgrid 
n = 100;
xaux = linspace(-pi,pi,n);
yaux =  linspace(0,2,n);
[x,y] = meshgrid(xaux,yaux);
nx = n-1; 
ny = n-1; 
h = 1/nx;

% generate right hand side vector b
b = zeros(n); 
b((-pi/2)<x & x<0.5*pi )=-cos(x((-pi/2)<x & x<0.5*pi ));
b(end,:) = [];
b(:,end) =[];
b = b(:);

% generate sparse matrix A 
I = speye(nx);
e = ones(nx,1);
T = spdiags([e -4*e circshift([2;e(2:end)],1)],[-1:1],nx,ny);
S = spdiags([e e e e],[-nx+1 -1 1 nx-1],nx,ny);
A = (kron(I,T) + kron(S,I))/h^2;

% Solve the linear system
u = zeros(n);
u_aux = A\b;
u(n:-1:2,1:n-1) = reshape(u_aux,ny,nx);
u(:,end)  = u(:,1); 



% plot the solution
figure(1);
clf; hold on; grid;
ma = max(max(u));
mi = min(min(u));
contourf(x,y,u,linspace(mi,ma,20));
title('numerical solution','fontsize',fsz);
xlabel('x','fontsize',fsz);
ylabel('y', 'fontsize',fsz);
