%% Direct numerical integration of the Saturn/Janus/Epimethius system

mu_s = 3.7931187e7; %km^3/s^2
mu_j = 0.1263; %km^3/s^2
mu_e = 0.0351; %km^3/s^2

%order: saturn, janus, epimetheus
x0 = [-1.865110336340190E+00, -9.865767926408326E+00, 2.457673802568250E-01;...
     -1.864109354923651E+00, -9.865703229168712E+00, 2.456363790938643E-01;...
     -1.865324283621047E+00, -9.864884006075078E+00, 2.453188536408405E-01].';

dx0 = [5.175097704745914E-03, -1.054146857424189E-03, -1.874911022818867E-04;...
      4.045755225099445E-03,  7.034449630331571E-03, -4.346218654110076E-03;...
      -3.713568443783068E-03, -2.500233487298903E-03, 1.399405483309747E-03].';

kmAU = 149597870.700;
mus = [mu_s;mu_j;mu_e];
mus = mus*86400^2/kmAU^3;

calcphi = @(z) acos(dot(z,[0,0,1])/norm(z));
calcth =  @(z) atan2(z(1)/norm(z(1:2)),-z(2)/norm(z(1:2)));
rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];


[t,x,dx] = nbodyVect(x0(:),dx0(:),mus,0:0.1:10*365.254,'c');
%[t,x,dx] = nbodyVect(x0(:),dx0(:),mus,0:0.001:100,'c');

xj = x(:,4:6) - x(:,1:3);
xe = x(:,7:9) - x(:,1:3);


figure(1)
clf()
hold on
jan = plot3(xj(1,1),xj(1,2),xj(1,3),'r.','MarkerSize',30);
eph = plot3(xe(1,1),xe(1,2),xe(1,3),'b.','MarkerSize',30);
jant = plot3(xj(1,1),xj(1,2),xj(1,3),'r-');
epht = plot3(xe(1,1),xe(1,2),xe(1,3),'b-');
view(3)
axis(reshape([min(xe);max(xe)]*1.1,6,1))

figure(2)
clf()
hold on
phi = calcphi(xj(1,:));
th = calcth(xj(1,:));
xjr = rotMatE2(-pi/2+phi)*rotMatE3(th)*xj(1,:).';
xer = rotMatE2(-pi/2+phi)*rotMatE3(th)*xe(1,:).';
janr = plot3(xjr(1),xjr(2),xjr(3),'r.','MarkerSize',30);
ephr = plot3(xer(1),xer(2),xer(3),'b.','MarkerSize',30);
axis([-1,1,-1,1]*0.0012)

for j = 2:length(t)
    set(jan,'XData',xj(j,1),'YData',xj(j,2),'ZData',xj(j,3));
    set(eph,'XData',xe(j,1),'YData',xe(j,2),'ZData',xe(j,3));
    set(jant,'XData',xj(1:j,1),'YData',xj(1:j,2),'ZData',xj(1:j,3));
    set(epht,'XData',xe(1:j,1),'YData',xe(1:j,2),'ZData',xe(1:j,3));
    
    phi = calcphi(xj(j,:));
    th = calcth(xj(j,:));
    xjr = rotMatE2(-pi/2+phi)*rotMatE3(-pi/2+th)*xj(j,:).';
    xer = rotMatE2(-pi/2+phi)*rotMatE3(-pi/2+th)*xe(j,:).';
    set(janr,'XData',xjr(1),'YData',xjr(2),'ZData',xjr(3));
    set(ephr,'XData',xer(1),'YData',xer(2),'ZData',xer(3));
    pause(1e-6);
end



