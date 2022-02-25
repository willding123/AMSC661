% Generate meshes 
% fd=@(p) ddiff(drectangle(p,0,3,0,3),dcircle(p,1.5,1.5,1));
th = 0:pi/100:2*pi;
r=1;
x = r * cos(th) + 1.5;
y = r * sin(th) + 1.5;
fixed_pt_c = [x;y]'; % fixed points on circle

fd = @(p) drectangle(p,0,3,0,3);
[p,t]=distmesh2d(fd,@huniform,0.1,[0,0;3,3],[fixed_pt_c; 0,0;0,3; 3,3; 3,0]); % generate mesh 

%% Boundary Condition
% Neumann BC 
% ind = zeros(lc1,1);
% for i = 1 : lc1
%     ind(i) = find(pts(:,1) == c1(i,1) & pts(:,2) == c1(i,2));
% end
% neumann = [ind,circshift(ind,[-1, 0])];
tol = 1e-2; 
up = find(p(:,2)>= 3-tol);% choose points in a small neighborhood of boundary
up = [up, circshift(up, [-1,0])];
up(length(up),:) =[];
dwn = find(p(:,2) <= tol);
dwn = [dwn, circshift(dwn, [-1,0])];
dwn(length(dwn),:) =[];
neumann = [up; dwn]; 

% ind = find(p(:,2)>=2.99); 
% neu = [ind,circshift(ind,[-1, 0])];
% neu(length(ind), :) = [];
% ind1 = find(p(:,2)<=0.01); 
% neu1 = [ind1,circshift(ind,[-1, 0])];
% neu1(length(ind1), :) = [];
% neumann = [neu; neu1];

% dirichlet
dirichlet = find(p(:,1)<=tol); 
dirichlet = [dirichlet; find(p(:,1)>=3-tol)];

% call FEM
u = FEM2D(p,t,neumann,dirichlet);

% calculate current density
Ntri = size(t,1); 
abs_current_centers = zeros(Ntri, 1);
for j = 1:Ntri % for all triangles    
  x = mean(p(t(j,:),1));
  y = mean(p(t(j,:),2));
  a = calc_a(x,y,.8,1); 
  G = Calc_G(p(t(j,:),:));
  abs_current_centers(j) =abs( sum(-a.*u(t(j,:))'*G(:,1)));
  % stima3_2 computes M = 0.5*|T_j|*G*G';
end

Npts = size(p,1);
abs_current_verts = zeros(Npts,1);
count_tri = zeros(Npts,1);
 for j = 1:Ntri
     abs_current_verts(t(j,:)) = abs_current_verts(t(j,:)) ...
           + abs_current_centers(j);
      count_tri(t(j,:)) = count_tri(t(j,:)) + 1;
 end
abs_current_verts = abs_current_verts./count_tri;
  
% graphic representation
figure(1);
trisurf(t,p(:,1),p(:,2),full(u),'facecolor','interp')
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
title('Voltage a1=0.8, a2=1')
set(gca,'Fontsize',20);
colorbar


figure(2);
trisurf(t,p(:,1),p(:,2),full(abs_current_verts),'facecolor','interp')
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
title("Current density a1=0.8, a2=1")
set(gca,'Fontsize',20);
colorbar



%% FEM
function u = FEM2D(pts,tri,neumann,dirichlet)
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,dirichlet); %mesh points with unknown values of u
A = sparse(Npts,Npts);
b = sparse(Npts,1);

%% Assembly
%% The Stiffness matrix
for j = 1:Ntri % for all triangles    
  x = mean(pts(tri(j,:),1));
  y = mean(pts(tri(j,:),2));
  a = calc_a(x,y,.8,1); 
  A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + a.*stima3(pts(tri(j,:),:));
  % stima3_2 computes M = 0.5*|T_j|*G*G';
end

%% The Right-hand side, i.e., the load vector
% Volume Forces
for j = 1:Ntri
  b(tri(j,:)) = 0;  % for the case where f = 0
end

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:)) = b(neumann(j,:)) + norm(pts(neumann(j,1),:)- ...
      pts(neumann(j,2),:)) * myg(sum(pts(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(Npts,1);
u(dirichlet) = myu_d(pts(dirichlet,:));
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


function G = Calc_G(verts)
G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
end

function a = calc_a(x,y, a1, a2)
if (x-1.5).^2+ (y-1.5).^2<= 1
    a = a1 ;
else 
    a = a2;
end
end

