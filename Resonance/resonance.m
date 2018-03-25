%% 2:1 resonance
ssdat = load('solarSystemData.mat');

%gravitational parameters AU^3/day^2
mus = ssdat.mus(10);
muj = ssdat.mus(5);
mu = mus+muj;
mua = mus+eps(mus);

e = 0.0489;
a = 5.204;

P = 2*pi * sqrt(a^3/mu);
n = sqrt(mu/a^3);
tpj = 0;

t = linspace(0,10*P,5000);

Pa = P/2;
aa = ((Pa/2/pi)^2*mua)^(1/3);
ea = 0.1;
tpa = Pa/2;
na = sqrt(mua/aa^3);

M = mod(n*(t-tpj),2*pi);
E = invKepler(M,e);
r = a.*(1-e.*cos(E)); %AU
r_js = [a*(cos(E) - e),a*sqrt(1 - e^2)*sin(E)];

Ma = mod(na*(t-tpa),2*pi);
Ea = invKepler(Ma,ea);
r = aa.*(1-ea.*cos(Ea)); %AU
r_as = [aa*(cos(Ea) - ea),aa*sqrt(1 - ea^2)*sin(Ea)];

dt = mean(diff(t))*200;
figure(2)
clf
hold on
axis([-6,6,-6,6])
s = plot(0,0,'g.','MarkerSize',60);
j = plot(r_js(1,1),r_js(1,2),'r.','MarkerSize',40);
jl = plot(r_js(1,1),r_js(1,2),'r--');
a = plot(r_as(1,1),r_as(1,2),'b.','MarkerSize',20);
al = plot(r_as(1,1),r_as(1,2),'b--');
grid on

for k = 2:length(t)
    set(j,'XData',r_js(k,1),'YData',r_js(k,2))
    set(a,'XData',r_as(k,1),'YData',r_as(k,2))
    set(jl,'XData',r_js(1:k,1),'YData',r_js(1:k,2))
    set(al,'XData',r_as(1:k,1),'YData',r_as(1:k,2))
    pause((t(k) - t(k-1))/dt)
end

%% badness
ssdat = load('solarSystemData.mat');

%gravitational parameters AU^3/day^2
mus = ssdat.mus(10);
muj = ssdat.mus(5);
mu = mus+muj;
mua = mus+eps(mus);

e = 0.0489;
a = 5.204;

P = 2*pi * sqrt(a^3/mu);
t = linspace(0,100*P,5000);

Pa = P/2;
aa = ((Pa/2/pi)^2*mua)^(1/3);
ea = 0.6;

%heliocentric planet position and velocity
r_js = [a*(1 - e),0,0];
v_js = [0,sqrt(mu*a)*sqrt(1 - e^2)*1./(a*(1-e)),0];
% 
% r_as = [aa*(1 - ea),0,0];
% v_as = [0,sqrt(mua*aa)*sqrt(1 - ea^2)*1./(aa*(1-ea)),0];
r_as = [aa*(-1 - ea),0,0];
v_as = [0,-sqrt(mua*aa)*sqrt(1 - ea^2)*1./(aa*(1+ea)),0];


%planet and sun vectors wrt barycenter
r_jb = mus*r_js/mu;
r_sb = -muj*r_js/mu;

v_jb = mus*v_js/mu;
v_sb = -muj*v_js/mu;

r_ab = r_as - r_sb;
v_ab = v_as - v_sb;

[T,Y,DY] = nbodyVect([r_jb,r_ab,r_sb].',[v_jb,v_ab,v_sb].',[muj,eps(mus),mus],t,'c');

dt = mean(diff(T))*1000;
figure(1)
clf
hold on
axis([min(Y(:)),max(Y(:)),min(Y(:)),max(Y(:))])
s = plot(Y(1,end-2),Y(1,end-1),'g.','MarkerSize',60);
j = plot(Y(1,1),Y(1,2),'r.','MarkerSize',40);
jl = plot(Y(1,1),Y(1,2),'r--');
a = plot(Y(1,4),Y(1,5),'b.','MarkerSize',20);
al = plot(Y(1,4),Y(1,5),'b--');

for k = 2:length(T)
    set(s,'XData',Y(k,end-2),'YData',Y(k,end-1))
    set(j,'XData',Y(k,1),'YData',Y(k,2))
    set(a,'XData',Y(k,4),'YData',Y(k,5))
    set(jl,'XData',Y(1:k,1),'YData',Y(1:k,2))
    set(al,'XData',Y(1:k,4),'YData',Y(1:k,5))
    pause((T(k) - T(k-1))/dt)
end


%% periapse advance
Es = linspace(0,2*pi,500);
%[r,v] = orbElem2vec(a,e,E,I,omega,Omega,mus)
[r1,v1] = orbElem2vec(1/(1 - 0.8),0.8,0,0,Es,0,mua);
r1 = r1.';

figure(3)
clf
hold on
s = plot(0,0,'g.','MarkerSize',60);
a = plot(r1(1,1),r1(1,2),'b.','MarkerSize',10);
al = plot(r1(1,1),r1(1,2),'b');
ao = plot(r1(1,1),r1(1,2),'r--');
axis([-10,10,-10,10])
axis square

for j = 1:10
for k = 1:length(Es)
    [rt,vt] = orbElem2vec(1/(1 - 0.8),0.8,Es,0,ones(size(Es))*Es(k),0,mua);
    set(a,'XData',r1(k,1),'YData',r1(k,2))
    set(al,'XData',r1(1:k,1),'YData',r1(1:k,2))
    set(ao,'XData',rt(1,:),'YData',rt(2,:))
    axis([-10,10,-10,10])
axis square
    pause(0.001);
end
end



