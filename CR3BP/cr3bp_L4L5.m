function [t,z] = cr3bp_L4L5(varargin)

% %large tadpole
%cr3bp_L4L5(0.001,[0.008,0,0.008,0],'animate',true)

% %approx horseshoe
%cr3bp_L4L5(0.000953875,[-1.475726125,0,-sqrt(3)/2,-0.06118],'tend',59*pi,'animate',true)

% %horseshoe
%cr3bp_L4L5(0.000953875,[-1.526496125,0,-sqrt(3)/2,0.04032],'tend',59*pi,'animate',true)


% input parsing
p = inputParser;
addOptional(p,'mu',0.001,@isnumeric);
addOptional(p,'z0offset',[0.0065,0,0.00065,0],@(x) isnumeric(x) && length(x)==4);
addOptional(p,'tend',32*pi,@isnumeric);
addParameter(p,'animate',false,@islogical);
addParameter(p,'fps',30,@islogical);
addOptional(p,'outfile','',@ischar);
parse(p,varargin{:});

% define system and calculate Lagrange points
mu = p.Results.mu;
z0 = p.Results.z0offset + [0.5-mu, 0, sqrt(3)/2, 0];

%co-linear L-point locs:
f = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
L3 = fsolve(f,-mu-0.1,optimoptions(@fsolve,'Display','none'));
L2 = fsolve(f,1-mu+0.1,optimoptions(@fsolve,'Display','none'));
L1 = fsolve(f,0.1,optimoptions(@fsolve,'Display','none'));

%energy of L-points
C4 = -1/2*(3+mu - mu^2);
Cl123 = @(x0) -(x0.^2)/2 - ((1-mu)./sqrt((x0+mu).^2) + mu./sqrt((x0 - (1-mu)).^2));


tstep = p.Results.tend/30/p.Results.fps;
[t,z] = ode113(@cr3bp_eom,0:tstep:p.Results.tend,z0,odeset('AbsTol',1e-16,'RelTol',1e-13));
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

szs = get(0, 'MonitorPositions');
sz = szs(end,:);

f = figure(1);
clf()
if ~isempty(p.Results.outfile)
    set(f,'Position',[sz(1), sz(2), sz(3)/2, sz(4)/2])
else
    set(f,'Position',sz)
end
subplot(1,2,1)
hold on
ax1 = gca;
set(ax1,'FontName','Times','FontSize',16)
contourf(ax1,X,Y,U,[Cl123(L1),Cl123(L2),Cl123(L3),C4])
plot(ax1,[-mu,1-mu],[0,0],'r.','MarkerSize',30)
plot(ax1,[0.5-mu,0.5-mu,L1,L2,L3],[sqrt(3)/2,-sqrt(3)/2,0,0,0],'k.','MarkerSize',20)
axis(ax1,'equal')
xlabel(ax1,'$\mathbf{\hat{e}}_r$','Interpreter','LaTex')
ylabel(ax1,'$\mathbf{\hat{e}}_\theta$','Interpreter','LaTex')

subplot(1,2,2)
hold on
ax2 = gca;
set(ax2,'FontName','Times','FontSize',16)
axis(ax2,'equal')
xlabel(ax2,'$\mathbf{\hat{e}}_1$','Interpreter','LaTex')
ylabel(ax2,'$\mathbf{\hat{e}}_2$','Interpreter','LaTex')

if p.Results.animate
    op = plot(ax1,x(1),y(1),'k','Linewidth',3);

    m1 = plot(ax2,xm1i(1),ym1i(1),'r.','MarkerSize',60);
    m2 = plot(ax2,xm2i(1),ym2i(1),'r.','MarkerSize',60);
    ppi = plot(ax2,xi(1),yi(1),'k.','MarkerSize',30);
    m1p = plot(xm1i(1),ym1i(1),'r');
    m2p = plot(xm2i(1),ym2i(1),'r','Linewidth',2);
    opi = plot(xi(1),yi(1),'k--');
    axis(ax2,[min([xm1i;xm2i;xi]),max([xm1i;xm2i;xi]),min([ym1i;ym2i;yi]),max([ym1i;ym2i;yi])])
    
    if ~isempty(p.Results.outfile)
        set(ax1,'nextplot','replacechildren');
        set(ax2,'nextplot','replacechildren');
        set(f,'Visible','off','Renderer','zbuffer')
        vidObj = VideoWriter(p.Results.outfile,'MPEG-4');
        vidObj.Quality = 100;
        vidObj.FrameRate = p.Results.fps;
        open(vidObj);
    end
    
    for j=2:length(x)
        set(op,'XData',x(1:j),'YData',y(1:j))
        set(opi,'XData',xi(1:j),'YData',yi(1:j))
        set(m1p,'XData',xm1i(1:j),'YData',ym1i(1:j))
        set(m2p,'XData',xm2i(1:j),'YData',ym2i(1:j))
        set(m1,'XData',xm1i(j),'YData',ym1i(j))
        set(m2,'XData',xm2i(j),'YData',ym2i(j))
        set(ppi,'XData',xi(j),'YData',yi(j))
        if ~isempty(p.Results.outfile)
            writeVideo(vidObj,getframe(f));
            disp(j/length(x)*100)
        else
            pause(1/p.Results.fps);
        end
    end
    if ~isempty(p.Results.outfile)
        close(vidObj);
        set(f,'Visible','on')
    end
else
    plot(ax1,x,y,'w','Linewidth',3)

    plot(ax2,xi,yi,'k')
    plot(ax2,xm1i,ym1i,'b')
    plot(ax2,xm2i,ym2i,'r')
end

figure(2)
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