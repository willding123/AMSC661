function basic_iterative_methods()
% solve -\Delta u = f on [0,1]^2, u = 0 on the boundary,
% u_exact = x*y*(1-x)*(1-y)

p1 = 4;
for p=p1:8
% plot error vs iteration only for the coarses mesh 17-by-17
    p_index = p-p1+1;
    n = 2^p + 1;
    n1 = n - 1;
    t = linspace(0,1,n);
    [x,y] = meshgrid(t,t);
    f = 2*(x.^2 - x + y.^2 - y);
    h(p_index) = 1/(n - 1);
    h2 = h(p_index)^2;
    u = zeros(n);
    u_exact = x.*y.*(1 - x).*(1 - y);

    itermax = 100;
    itermin = 50;
    ii = itermin : itermax;
    om = 1.5:0.01:2;   
for o =1:length(om)
    % SOR
    u = zeros(n);
%     om = 2/(1 + sin(pi*h(p_index))) + eps;
    for k = 1 : itermax
    for i = 2 : n1
        for j = 2 : n1
            u(i,j) = (1 - om(o))*u(i,j) + om(o)*0.25*(u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - h2*f(i,j));
        end
    end
    esor(k) = max(max(abs(u - u_exact)));
    end
    z = log(esor(itermin : itermax));
    pol = polyfit(ii,z,1);
    rsor(o) = -pol(1);
end
%% plot convergence rate
figure; clf; grid;
hold on;
plot(om,rsor','k','Linewidth',2);
set(gca,'Fontsize',20);
xlabel('w');
ylabel('Convergence Rate');
title('N=2^'+string(p) +"+1")

end
ome = [1.68,1.83, 1.91, 1.97, 1.99];
figure(1); clf; grid;
hold on;
plot(h,ome,'k','Linewidth',2);
hold on 
plot(h, 2./(1 + sin(pi.*h)) + eps)
set(gca,'Fontsize',20);
xlabel('h');
ylabel('w');
legend('Experimental', 'Theoretical')
% title('N=2^'+string(p) +"+1")



end


