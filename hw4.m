%% Generate mesh 

% N=[18 36 72];
% error = []; 
% hval = [];
% for k = 1:length(N)
N =72;
h = pi/N;
hval(k) = h;
th = 0:h:2*pi;
r=1;
x = r * cos(th) ;
y = r * sin(th) ;
fixed_pt_c = [x;y]'; % fixed points on circle
fd = @(p) ddiff(dcircle(p,0,0,2), dcircle(p,0,0,1));

% generate equidistant points on a circle 
r1 =2; 
a = 2*asin(h/(2*r1));
x1 =[];
y1= []; 
for i=1:floor(2*pi/a)
    x1 = [x1; r1*cos(i*a)];
    y1 = [y1; r1*sin(i*a)];
end 
fixed_pt_c1 =[x1 y1]; 
[p,t]=distmesh2d(fd,@huniform,h,[-2,-2;2,2],[fixed_pt_c; fixed_pt_c1 ]); % generate mesh 

% Dirichlet 

tol =1e-4; 
inner = find(p(:,1 ).^2+ p(:,2 ).^2<=(r+tol)^2 & p(:,1 ).^2+ p(:,2 ).^2>=(r-tol)^2); % innner circle
outer = find(p(:,1 ).^2+ p(:,2 ).^2<=(r1+tol)^2 & p(:,1 ).^2+ p(:,2 ).^2>=(r1-tol)^2); % outer circle 
dirichlet = [inner; outer];
neumann = []; 
u = FEM2D(p,t,neumann,dirichlet);

% graphic representation
figure();

hold on
axis ij
view(2)
trisurf(t,p(:,1),p(:,2),full(u),'facecolor','interp')

xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
colorbar

% plot exact solution against computed
figure();
u_exact = @(r)-0.25*r.^2 + 0.75*log(r)/log(2) + 0.25;
R = sqrt(p(:,1).^2 + p(:,2).^2);
[rsort,isort] = sort(R,'ascend');
usort = u(isort);
ue = u_exact(rsort);
plot(rsort,usort, '--ob','Linewidth',2);
hold on
plot(rsort,ue,'Linewidth',2);
xlabel('r','Fontsize',20);
ylabel('u','Fontsize',20);
legend("computed solution","exact solution")
set(gca,'Fontsize',20);

% error(k) =max(abs(usort-ue));

% end 
% figure; hold on;
% plot(hval,error,'.','Markersize',20);
% x = log(hval);
% y = log(error);
% p = polyfit(x,y,1);
% pval = polyval(p,x);
% plot(hval,exp(pval),'Linewidth',2);
% set(gca,'YScale','log','XScale','log','Fontsize',20);
% xlabel('h','Fontsize',20);
% ylabel('max abs error','Fontsize',20);

%% 
function u = FEM2D(pts,tri,neumann,dirichlet)
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,dirichlet); %mesh points with unknown values of u
A = sparse(Npts,Npts);
b = sparse(Npts,1);

%% Assembly
%% The Stiffness matrix
for j = 1:Ntri % for all triangles    
  A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + stima3(pts(tri(j,:),:));
  % stima3_2 computes M = 0.5*|T_j|*G*G';
end

%% The Right-hand side, i.e., the load vector
% Volume Forces
for j = 1:Ntri
  b(tri(j,:)) = b(tri(j,:)) + ...
det([1 1 1; pts(tri(j,:),:)'])/6;  
end

% f(sum(coordinates(elements3(j,:),

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:)) = b(neumann(j,:)) + norm(pts(neumann(j,1),:)- ...
      pts(neumann(j,2),:)) * myg(sum(pts(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(Npts,1);
u(dirichlet) = 0;
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

end

%%
function DirichletBoundaryValue = myu_d(x)
xmin = min(x(:,1));
xmax = max(x(:,1));
midx = 0.5*(xmin + xmax);
DirichletBoundaryValue =  0.5 * (sign(x(:,1) - midx) + 1);
end

%%
function Stress = myg(x)
Stress = zeros(size(x,1),1);
end

%
function M = stima3(verts)
G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
M = 0.5*det([ones(1,3);verts']) * G * G';
end

