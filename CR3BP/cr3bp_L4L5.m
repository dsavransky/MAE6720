function [t,z,pmap] = cr3bp_L4L5(varargin)

% %large tadpole
%cr3bp_L4L5(0.001,[0.008,0,0.008,0],'animate',true)

% %approx horseshoe
%cr3bp_L4L5(0.000953875,[-1.475726125,0,-sqrt(3)/2,-0.06118],200,'animate',true)

% %horseshoe
%cr3bp_L4L5(0.000953875,[-1.526496125,0,-sqrt(3)/2,0.04032],200,'animate',true)


% input parsing
p = inputParser;
addOptional(p,'mu',0.001,@isnumeric);
addOptional(p,'z0offset',[0.0065,0,0.00065,0],@(x) isnumeric(x) && length(x)==4);
addOptional(p,'tend',100,@isnumeric);
addParameter(p,'animate',false,@islogical);
parse(p,varargin{:});

% define system and calculate Lagrange points
mu = p.Results.mu;
z0 = p.Results.z0offset + [0.5-mu, 0, sqrt(3)/2, 0];

%co-linear L-point locs:
f = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
L3 = fsolve(f,-mu-0.1);
L2 = fsolve(f,1-mu+0.1);
L1 = fsolve(f,0.1);

%energy of L-points
C4 = -1/2*(3+mu - mu^2);
Cl123 = @(x0) -(x0.^2)/2 - ((1-mu)./sqrt((x0+mu).^2) + mu./sqrt((x0 - (1-mu)).^2));


[t,z] = ode113(@cr3bp_eom,[0:0.01:p.Results.tend],z0,odeset('AbsTol',1e-16,'RelTol',1e-13));
x = z(:,1);
xd = z(:,2);
y = z(:,3);
yd = z(:,4);
r1 = sqrt((x+mu).^2+y.^2);
r2 = sqrt((x - (1-mu)).^2+y.^2);
U = -(x.^2+y.^2)/2 - ((1-mu)./r1 + mu./r2);
C = (xd.^2 + yd.^2)/2. + U;

tmp = linspace(-1.5,1.5,500);
[X,Y] = meshgrid(tmp,tmp);
r1 = sqrt((X+mu).^2+Y.^2);
r2 = sqrt((X - (1-mu)).^2+Y.^2);
U = -(X.^2+Y.^2)/2 - ((1-mu)./r1 + mu./r2);

%inertial frame
xi = x.*cos(t)-y.*sin(t);
yi = x.*sin(t) + y.*cos(t);

xm1i = -mu*cos(t);
ym1i = -mu*sin(t);
xm2i = (1-mu)*cos(t);
ym2i = (1-mu)*sin(t);


figure(1)
clf()
hold on
set(gca,'FontName','Times','FontSize',16)
contourf(X,Y,U,[Cl123(L1),Cl123(L2),Cl123(L3),C4])
plot([-mu,1-mu],[0,0],'r.','MarkerSize',30)
plot([0.5-mu,0.5-mu,L1,L2,L3],[sqrt(3)/2,-sqrt(3)/2,0,0,0],'k.','MarkerSize',20)
figure(2)
clf()
hold on
set(gca,'FontName','Times','FontSize',16)
if p.Results.animate
    figure(1)
    op = plot(x(1),y(1),'k','Linewidth',3);
    figure(2)
    m1 = plot(xm1i(1),ym1i(1),'r.','MarkerSize',60);
    m2 = plot(xm2i(1),ym2i(1),'r.','MarkerSize',50);
    ppi = plot(xi(1),yi(1),'k.','MarkerSize',30);
%     m1p = plot(xm1i(1),ym1i(1),'r');
%     m2p = plot(xm2i(1),ym2i(1),'r');
%     opi = plot(xi(1),yi(1),'k');
    axis([min([xm1i;xm2i;xi]),max([xm1i;xm2i;xi]),min([ym1i;ym2i;yi]),max([ym1i;ym2i;yi])])
    for j=2:length(x)
        set(op,'XData',x(1:j),'YData',y(1:j))
%        set(opi,'XData',xi(max([1,j-100]):j),'YData',yi(max([1,j-100]):j))
%         set(m1p,'XData',xm1i(1:j),'YData',ym1i(1:j))
%         set(m2p,'XData',xm2i(1:j),'YData',ym2i(1:j))
        set(m1,'XData',xm1i(j),'YData',ym1i(j))
        set(m2,'XData',xm2i(j),'YData',ym2i(j))
        set(ppi,'XData',xi(j),'YData',yi(j))
        pause((t(j) - t(j-1))/10);
    end
else
    figure(1)
    plot(x,y,'k')
    figure(2)
    plot(xi,yi,'k')
    plot(xm1i,ym1i,'b')
    plot(xm2i,ym2i,'r')
end

figure(3)
clf()
plot(t,C)


    function dz = cr3bp_eom(t,z)
        
        x = z(1); xd = z(2); y = z(3); yd = z(4);
        
        xdd = x-mu.*1.0./((mu+x-1.0).^2+y.^2).^(3.0./2.0).*(mu.*2.0+x.*2.0-2.0).*...
            (1.0./2.0)+(mu.*2.0+x.*2.0).*(mu-1.0).*1.0./...
            ((mu+x).^2+y.^2).^(3.0./2.0).*(1.0./2.0) + 2*yd;
        ydd = y-mu.*y.*1.0./((mu+x-1.0).^2+y.^2).^(3.0./2.0)+y.*(mu-1.0).*1.0./...
            ((mu+x).^2+y.^2).^(3.0./2.0) - 2*xd;
        
        dz = [xd;xdd;yd;ydd];
        
    end


end