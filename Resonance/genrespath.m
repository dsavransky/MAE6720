function genrespath(ea,f)
% genrespath(ea, f) animates the paths of two bodies in resonant motion
% with the outer body on an orbit equivalent to Jupiter's, and the inner 
% body on an orbit with eccentricity ea, and period ratio to jupiter of f.

% Eaxmples:
%   genrespath(0.3,2/1)
%   genrespath(0.3,4/3)
%   genrespath(0.3,7/5)
%   genrespath(0.3,3/4)

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

fps = 60;
nframes = 10*fps;

t = linspace(0,5*P,nframes);

Pa = P/f;
aa = ((Pa/2/pi)^2*mua)^(1/3);
tpa = 0;
na = sqrt(mua/aa^3);

M = mod(n*(t-tpj),2*pi);
E = invKepler(M,e);
r_js = [a*(cos(E) - e),a*sqrt(1 - e^2)*sin(E)];

Ma = mod(na*(t-tpa),2*pi);
Ea = invKepler(Ma,ea);
r_as = [aa*(cos(Ea) - ea),aa*sqrt(1 - ea^2)*sin(Ea)];

rc = zeros(size(r_as));
for k = 1:length(t)
    rc(k,:) = ([cos(n*t(k)) sin(n*t(k)); -sin(n*t(k)) cos(n*t(k))]*r_as(k,:).').'; 
end

dt = mean(diff(t))*500;
figure(2)
clf
%set(2,'Position',[200,200,1120,420])
subplot(1,2,1)
hold on
tmp = max(abs([r_js(:);r_as(:)]));
axis equal
axis([-tmp,tmp,-tmp,tmp])
plot(0,0,'g.','MarkerSize',60);
j = plot(r_js(1,1),r_js(1,2),'r.','MarkerSize',40);
jl = plot(r_js(1,1),r_js(1,2),'r--');
a = plot(r_as(1,1),r_as(1,2),'b.','MarkerSize',20);
al = plot(r_as(1,1),r_as(1,2),'b--');
grid on


subplot(1,2,2)
hold on
tmp = max(abs(rc(:)));
axis equal
axis([-tmp,tmp,-tmp,tmp])
plot(0,0,'g.','MarkerSize',60);
ac = plot(rc(1,1),rc(1,2),'b.','MarkerSize',20);
acl = plot(rc(1,1),rc(1,2),'b--');

for k = 2:length(t)
    set(j,'XData',r_js(k,1),'YData',r_js(k,2))
    set(a,'XData',r_as(k,1),'YData',r_as(k,2))
    set(ac,'XData',rc(k,1),'YData',rc(k,2))
    set(jl,'XData',r_js(1:k,1),'YData',r_js(1:k,2))
    set(al,'XData',r_as(1:k,1),'YData',r_as(1:k,2))
    set(acl,'XData',rc(1:k,1),'YData',rc(1:k,2))
    pause(1/fps);
    %pause((t(k) - t(k-1))/dt*2)
end

