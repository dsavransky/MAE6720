% Animations of the various natural motions in the Euler-Hill frame

%% Case 1
n = 1;
r1 = 1;
x0 = 0.01;
y0 = 0;
yd0 = -3*n*x0/2;
t = linspace(0,6*pi,1000);
x = x0*ones(size(t));
y = -3*n*t*x0/2 + y0;
M = mod(n*t,2*pi);
reforb = [r1*cos(M);r1*sin(M)];
hillorb_rot = [x;y];
hillorb = hillorb_rot;
for j = 1:length(t)
    hillorb(:,j) = [cos(M(j)), -sin(M(j)); sin(M(j)) ,cos(M(j))]*hillorb_rot(:,j);
end
hillorb = hillorb + reforb;

figure(11)
clf
set(gca,'Xlim',[-1.25,1.25],'YLim',[-1.25,1.25])
axis equal
hold on
plot(0,0,'k.','MarkerSize',50);
msr = plot(reforb(1,1),reforb(2,1),'b.','MarkerSize',20);
ms = plot(hillorb(1,1),hillorb(2,1),'r.','MarkerSize',40);
ps = plot(reforb(1,1),reforb(2,1),'b--',hillorb(1,1),hillorb(2,1),'r--');
hold off
set(gca,'Xlim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],...
    'YLim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],'YTick',[],'Xtick',[])
for j = 2:length(t)
    set(ms,'XData',hillorb(1,j),'YData',hillorb(2,j));
    set(msr,'XData',reforb(1,j),'YData',reforb(2,j));
    set(ps(1),'XData',reforb(1,1:j),'YData',reforb(2,1:j));
    set(ps(2),'XData',hillorb(1,1:j),'YData',hillorb(2,1:j));
    pause(0.01);
end

%% case 2
n = 1;
r1 = 1;
x0 = 0.05;
y0 = 0.05;
yd0 = -2*n*x0;
t = linspace(0,6*pi,1000);
x = x0*cos(n*t);
y = -2*x0*sin(n*t) + y0;
M = n*t;
reforb = [r1*cos(M);r1*sin(M)];
hillorb_rot = [x;y];
hillorb = hillorb_rot;
for j = 1:length(t)
    hillorb(:,j) = [cos(M(j)), -sin(M(j)); sin(M(j)) ,cos(M(j))]*hillorb_rot(:,j);
end
hillorb = hillorb + reforb;
figure(12)
plot(x,y)
axis equal

figure(13)
clf
set(gca,'Xlim',[-1.25,1.25],'YLim',[-1.25,1.25])
axis equal
hold on
plot(0,0,'k.','MarkerSize',50);
msr = plot(reforb(1,1),reforb(2,1),'b.','MarkerSize',20);
ms = plot(hillorb(1,1),hillorb(2,1),'r.','MarkerSize',40);
ps = plot(reforb(1,1),reforb(2,1),'b--',hillorb(1,1),hillorb(2,1),'r--');
hold off
set(gca,'Xlim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],...
    'YLim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],'YTick',[],'Xtick',[])
for j = 2:length(t)
    set(ms,'XData',hillorb(1,j),'YData',hillorb(2,j));
    set(msr,'XData',reforb(1,j),'YData',reforb(2,j));
    set(ps(1),'XData',reforb(1,1:j),'YData',reforb(2,1:j));
    set(ps(2),'XData',hillorb(1,1:j),'YData',hillorb(2,1:j));
    pause(0.01);
end

%% case 3 - drift - x0=y0=xd0=0, yd0 >0
n = 3;
r1 = 1;
x0 = 0;
y0 = 0;
yd0 = 0.01;
t = linspace(0,6*pi,1000);
x = 2*(1-cos(n*t))*yd0/n;
y = yd0/n*(4*sin(n*t) - 3*n*t);
M = n*t;
reforb = [r1*cos(M);r1*sin(M)];
hillorb = reforb+[x;y];

figure(12)
plot(x,y)
axis equal

figure(13)
clf
set(gca,'Xlim',[-1.25,1.25],'YLim',[-1.25,1.25])
axis equal
hold on
plot(0,0,'k.','MarkerSize',30);
ms = plot(hillorb(1,1),hillorb(2,1),'r.','MarkerSize',20);
ps = plot(reforb(1,1),reforb(2,1),'b--',hillorb(1,1),hillorb(2,1),'r--');
hold off
set(gca,'Xlim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],'YLim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))])
for j = 2:length(t)
    set(ms,'XData',hillorb(1,j),'YData',hillorb(2,j));
    set(ps(1),'XData',reforb(1,1:j),'YData',reforb(2,1:j));
    set(ps(2),'XData',hillorb(1,1:j),'YData',hillorb(2,1:j));
    pause(0.001);
end

%% Gas drag 1

n = 1;
r1 = 1;
x0 = 0;
y0 = 0;
D = 0.005;
t = linspace(0,2*pi,100);
x = -2/n*D*t;
y = 3*D*t.^2/2;
M = n*t;
reforb = [r1*cos(M);r1*sin(M)];
hillorb = reforb+[x;y];

figure(12)
plot(x,y)
axis equal

figure(13)
clf
set(gca,'Xlim',[-1.25,1.25],'YLim',[-1.25,1.25])
axis equal
hold on
plot(0,0,'k.','MarkerSize',30);
ms = plot(hillorb(1,1),hillorb(2,1),'r.','MarkerSize',20);
ps = plot(reforb(1,1),reforb(2,1),'b--',hillorb(1,1),hillorb(2,1),'r--');
hold off
set(gca,'Xlim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],'YLim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))])
for j = 2:length(t)
    set(ms,'XData',hillorb(1,j),'YData',hillorb(2,j));
    set(ps(1),'XData',reforb(1,1:j),'YData',reforb(2,1:j));
    set(ps(2),'XData',hillorb(1,1:j),'YData',hillorb(2,1:j));
    pause(0.1);
end