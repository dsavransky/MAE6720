function [t,z,pmap] = cr3bp_poincare(varargin)

% %stable periodic
% [t,z,pmap] = cr3bp_poincare('animate',false,'X0',0.0935,'C0',-5,'n_poincare',5000);
% [t,z,pmap] = cr3bp_poincare('animate',true,'X0',0.0935,'C0',-5,'tend',2*pi,'outfile','cr3b_poincare_x0_0p0935_C0_n5.mp4');


% %stable resonant
% [t,z,pmap] = cr3bp_poincare('animate',false,'X0',0.1,'n_poincare',5000)
% [t,z,pmap] = cr3bp_poincare('animate',true,'X0',0.1,'tend',12*pi,'fpo',300,'outfile','cr3b_poincare_x0_0p1_C0_L2.mp4')

% %chaotic
% [t,z,pmap] = cr3bp_poincare('animate',false,'X0',0.18,'n_poincare',2500)
% [t,z,pmap] = cr3bp_poincare('animate',true,'X0',0.18,'tend',12*pi,'fpo',150,'outfile','cr3b_poincare_x0_0p18_C0_L2.mp4')

%stable quasi-periodic
% [t,z,pmap] = cr3bp_poincare('animate',false,'X0',0.9,'n_poincare',5000)
% [t,z,pmap] = cr3bp_poincare('animate',true,'X0',0.9,'tend',12*pi)


% [t,z,pmap] = cr3bp_poincare('animate',false,'n_poincare',2500)

%define system and calculate Lagrange points
mu = 1/82.3; %earth-moon
%co-linear L-point locs:
f = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
L3 = fsolve(f,-mu-0.1,optimoptions(@fsolve,'Display','none'));
L2 = fsolve(f,1-mu+0.1,optimoptions(@fsolve,'Display','none'));
L1 = fsolve(f,0.1,optimoptions(@fsolve,'Display','none'));

%energy of L-points
C4 = -1/2*(3- mu + mu^2);
Cl123 = @(x0) -(x0.^2)/2 - ((1-mu)./sqrt((x0+mu).^2) + mu./sqrt((x0 - (1-mu)).^2));

%input parsing
p = inputParser;
addParameter(p,'n_poincare',0,@(x) isnumeric(x) && x>=0);
addParameter(p,'x0',L2,@isnumeric);
addParameter(p,'C0',Cl123(L2),@isnumeric);
addParameter(p,'animate',false,@islogical);
addOptional(p,'tend',50,@isnumeric);
addOptional(p,'fpo',30,@isnumeric);
addOptional(p,'outfile','',@ischar);

parse(p,varargin{:});
n_poincare = p.Results.n_poincare;
x0 = p.Results.x0;
C0 = p.Results.C0;

%orbit initial condition
y0 = 0;
xd0 = 0;
r1 = sqrt((x0+mu).^2+y0.^2);
r2 = sqrt((x0 - (1-mu)).^2+y0.^2);
U = -(x0.^2+y0.^2)/2 - ((1-mu)./r1 + mu./r2);
yd0 = sqrt(2*(C0-U) - xd0.^2);

%let's figure out the orbital period if this were a two-body orbit
if x0 <= L1
    % 'orbiting' m1
    mup = 1 - mu;
    r = abs(x0+mu);
    v = yd0 + x0 + mu;
else
    % 'orbiting' m1
    mup = mu;
    r = abs(1 - mu - x0);
    v = yd0 + x0 - (1 - mu);
end
En = v^2/2 - mup/r;
a = -mup/2/En;
Tp = 2*pi/sqrt(mup)*a^(3/2);
%have one orbit take 1 second
tstep = Tp/p.Results.fpo;

z0 = [x0;xd0;y0;yd0];

if n_poincare > 0
    ze = z0.';
    te = 0;
    n = n_poincare;
    pmap = zeros(n,4);
    counter = 1;
    while counter <= n
        [~,~,te,ze] = ode113(@cr3bp_eom,[te(end),te(end)+100],ze(end,:),odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@cr3bp_ycross));
        pmap(counter,:) = ze(end,:);
        counter = counter + 1;
        if mod(counter,100) == 0,disp(counter); end
    end
    disp(te)
else
    pmap = [];
end

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


%sz = get(0,'Screensize');
szs = get(0, 'MonitorPositions');
sz = szs(end,:);

%trajectory plots
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
set(gca,'FontName','Times','FontSize',16)
contourf(ax1,X,Y,U,[Cl123(L1),Cl123(L2),Cl123(L3),C4])
plot(ax1,[-mu,1-mu],[0,0],'r.','MarkerSize',30)
plot(ax1,[0.5-mu,0.5-mu,L1,L2,L3],[sqrt(3)/2,-sqrt(3)/2,0,0,0],'k.','MarkerSize',20)
axis(ax1,'equal')
xlabel(ax1,'$\mathbf{\hat{e}}_r$','Interpreter','LaTex')
ylabel(ax1,'$\mathbf{\hat{e}}_\theta$','Interpreter','LaTex')

subplot(1,2,2)
hold on
ax2 = gca;
set(gca,'FontName','Times','FontSize',16)
axis(ax2,'equal')
xlabel(ax2,'$\mathbf{\hat{e}}_1$','Interpreter','LaTex')
ylabel(ax2,'$\mathbf{\hat{e}}_2$','Interpreter','LaTex')

if p.Results.animate
    op = plot(ax1,x(1),y(1),'k');

    m1 = plot(ax2,xm1i(1),ym1i(1),'r.','MarkerSize',40);
    m2 = plot(ax2,xm2i(1),ym2i(1),'r.','MarkerSize',40);
    m1p = plot(ax2,xm1i(1),ym1i(1),'r','Linewidth',2);
    m2p = plot(ax2,xm2i(1),ym2i(1),'r','Linewidth',2);
    opi = plot(ax2,xi(1),yi(1),'k');
    ppi = plot(ax2,xi(1),yi(1),'b.','MarkerSize',30);
    axis(ax2,[min([xm1i;xm2i;xi]),max([xm1i;xm2i;xi]),min([ym1i;ym2i;yi]),max([ym1i;ym2i;yi])])
    
    if ~isempty(p.Results.outfile)
        set(ax1,'nextplot','replacechildren');
        set(ax2,'nextplot','replacechildren');
        set(f,'Visible','off','Renderer','zbuffer')
        vidObj = VideoWriter(p.Results.outfile,'MPEG-4');
        vidObj.Quality = 100;
        vidObj.FrameRate = p.Results.fpo;
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
            pause(1/p.Results.fpo);
        end
    end
    if ~isempty(p.Results.outfile)
        close(vidObj);
        set(f,'Visible','on')
    end
else
    plot(ax1,x,y,'k')
    axis(ax2,[min([xm1i;xm2i;xi]),max([xm1i;xm2i;xi]),min([ym1i;ym2i;yi]),max([ym1i;ym2i;yi])])
    plot(ax2,xi,yi,'k')
    plot(ax2,xm1i,ym1i,'b')
    plot(ax2,xm2i,ym2i,'r')
end

% plot Jacobi constant throughout integration
figure(2)
clf()
plot(t,C)

% Poincare map
if n_poincare > 0
    figure(3)
    clf()
    plot(pmap(pmap(:,4)<=0,1),pmap(pmap(:,4)<=0,2),'.')
    set(gca,'FontName','Times','FontSize',16)
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\dot x$','Interpreter','Latex')
    title(sprintf('$x_0 = %1.4g, C = %1.4g$',x0,C0),'Interpreter','Latex')
end

    function dz = cr3bp_eom(~,z)
        
        x = z(1); xd = z(2); y = z(3); yd = z(4);
        
        xdd = x-mu.*1.0./((mu+x-1.0).^2+y.^2).^(3.0./2.0).*(mu.*2.0+x.*2.0-2.0).*...
            (1.0./2.0)+(mu.*2.0+x.*2.0).*(mu-1.0).*1.0./...
            ((mu+x).^2+y.^2).^(3.0./2.0).*(1.0./2.0) + 2*yd;
        ydd = y-mu.*y.*1.0./((mu+x-1.0).^2+y.^2).^(3.0./2.0)+y.*(mu-1.0).*1.0./...
            ((mu+x).^2+y.^2).^(3.0./2.0) - 2*xd;
        
        dz = [xd;xdd;yd;ydd];
        
    end

    function [value,isterminal,direction] = cr3bp_ycross(~,z)
        %find y = 0 crossings
        value = z(3);
        isterminal = 1;
        direction = 0;
    end

end