function myVcycle()
% solve -\nabla * [a(x,y) * \nabla u] = f
% (x,y) \in [0,1]^2; u = 0 on the boundary
% uexact = x.*(1 - x).*y.*(1 - y);
p = 10;
n = 2^p + 1;
h = 2^(-p);
n1 = n - 1;
% nsmooth = 1: V-cycle
% nsmooth = 2: W-cycle
nsmooth = 2;

tol = 1e-9;
Iinner = 2 : n1; % indices of the inner mesh points

t = linspace(0,1,n);
[x, y] = meshgrid(t,t);
uexact = x.*(1 - x).*y.*(1 - y);
a = 1 + x + 2*y.^2;
f = zeros(n);
f(Iinner,Iinner) = diff_oper(h,uexact,a);

Vcounter = 1; Vmax =40; 
u = zeros(n); % the initial guess for the solution
res(Vcounter) = max(max(abs(diff_oper(h,u,a) - f(Iinner,Iinner))));
err(Vcounter) = max(max(abs(u - uexact)));
fprintf('Vcounter = %d, res = %d, err = %d\n',Vcounter - 1,res(Vcounter),err(Vcounter));
tic
while res(Vcounter) > tol & Vcounter< Vmax
    u = Vcycle(u,a,f,p,nsmooth);
    Vcounter = Vcounter + 1;
    res(Vcounter) = max(max(abs(diff_oper(h,u,a) - f(Iinner,Iinner))));
    err(Vcounter) = max(max(abs(u - uexact)));
    fprintf('Vcounter = %d, res = %d, err = %d\n',Vcounter - 1,res(Vcounter),err(Vcounter));
end
CPUtime = toc;
fprintf('Multigrid V-cycle: CPUtime = %d\n',CPUtime);
figure;
clf;
hold on;
plot(0:Vcounter-1,res,'Linewidth',1,'marker','*','color','b');
set(gca,'YScale','log','Fontsize',24)
title('  Jacobi ','FontSize',24)
xlabel('V cycle #','Fontsize',20);
ylabel('Residual','FontSize',20);
grid
hold off



end

%%
function u = Vcycle(u,a,f,p,nsmooth)
N = 2^p + 1;
h = 2^(-p);
pcoarse = p - 1;
% smoothing
u = smooth(nsmooth,h,u,a,f);
r = zeros(N);
Iinner = 2 : N - 1;
r(Iinner,Iinner) = f(Iinner,Iinner) - diff_oper(h,u,a);
% restriction
Icoarse = 1 : 2 : N;
e = zeros(2^pcoarse + 1);
rcoarse = r(Icoarse,Icoarse);
acoarse = a(Icoarse,Icoarse);
H = 2*h;
if pcoarse == 1 % the 3-by-3 mesh is reached
    e(2,2) =  rcoarse(2,2)*2*H^2/(acoarse(1,2) + acoarse(2,1) + ...
        4*acoarse(2,2) + acoarse(2,3) + acoarse(3,2));
else
    e = Vcycle(e,acoarse,rcoarse,pcoarse,nsmooth);
end
% prolongation to the fine mesh
Ifine = 1 : N;
[Icmesh,Jcmesh] = meshgrid(Icoarse);
[Ifmesh,Jfmesh] = meshgrid(Ifine);
efine = interp2(Icmesh,Jcmesh,e,Ifmesh,Jfmesh);
% correct the solution
u = u + efine;
%postsmoothing
u = smooth(nsmooth,h,u,a,f);
end


%%
function Lu = diff_oper(h,u,a)
n = length(u);

    as = 0.5*(a + circshift(a,[1,0]));
    an = 0.5*(a + circshift(a,[-1,0]));
    aw = 0.5*(a + circshift(a,[0,1]));
    ae = 0.5*(a + circshift(a,[0,-1]));
    ap = aw + ae + as + an;
    I = 2 : n - 1;

Lu = (ap(I,I).*u(I,I) - as(I,I).*u(I - 1,I) ...
    - an(I,I).*u(I + 1,I) - aw(I,I).*u(I,I - 1)...
    - ae(I,I).*u(I,I + 1))/h^2;

end

%%

%%
function u = smooth(nsmooth,h,u,a,f)
n = length(u);
n1 = n - 1;
% coefficients for the finite difference scheme
% ((au_x)_x + (au_y)_y + f)(i,j) = 
% (1/(h^2))(aw*u(i-1,j) + ae*u(i+1,i) + as*u(i,j-1)+an*u(i,j+1) -
%            ap*u(i,j)) + f(i,j) = 0
as = 0.5*(a + circshift(a,[1,0]));
an = 0.5*(a + circshift(a,[-1,0]));
aw = 0.5*(a + circshift(a,[0,1]));
ae = 0.5*(a + circshift(a,[0,-1]));
ap = aw + ae + as + an;

h2 = h^2;
% % Gauss-Seidel
% for k = 1 : nsmooth
%     for i = 2 : n1
%         for j = 2 : n1
%             u(i,j) = (f(i,j)*h2 + u(i - 1,j)*as(i,j) + u(i + 1,j)*an(i,j)+...
%                 u(i,j - 1)*aw(i,j) + u(i,j + 1)*ae(i,j))/ap(i,j);
%         end
%     end
% end

% %weighted Jacobi 
% om =1; 
% Ip = [3 : n];
% Im = [1 : n - 2];
% I = 2:n1; 
% for k =1:nsmooth
% u(I,I) = (1 - om)*u(I,I) + om*(u(I,Ip).*ae(I,I) +u(I,Im).*aw(I,I) + u(Ip,2:n1).*as(I,I) + u(Im,I).*an(I,I) + h2*f(I,I))./(ap(I,I));
% end

% red black GS
for k =1:nsmooth
    for i = 2:n1
        for j = 2:n1
            if mod(i+j,2) == 0
                u(i,j) = (f(i,j)*h2 + u(i - 1,j)*as(i,j) + u(i + 1,j)*an(i,j)+...
                 u(i,j - 1)*aw(i,j) + u(i,j + 1)*ae(i,j))/ap(i,j); 
            end
        end
    end

    for i = 2:n1
        for j = 2:n1
            if mod(i+j,2) == 1
                u(i,j) = (f(i,j)*h2 + u(i - 1,j)*as(i,j) + u(i + 1,j)*an(i,j)+...
                 u(i,j - 1)*aw(i,j) + u(i,j + 1)*ae(i,j))/ap(i,j);
            end
        end
    end
end


end
