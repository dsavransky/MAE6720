% A comparison of the secular two-body solution with direct numerical 
% integration for Jupiter and Saturn

%% Numerical Integration of Jupiter/Saturn/Sun

mAU = 1.495978707e11; %m

% 1/1/1983
rsun = [5.835670348697200E-03; 7.803816747876102E-03; -2.393595507252648E-04];
vsun = [-6.753425418733883E-06; 5.280474347005359E-06; 1.223761876644164E-07];
musun = 1.3271244004193938e11; % km^3/s^2
musun = musun*(1000/mAU)^3 * (86400)^2;

rj = [-3.096932989848971E+00; -4.400340838525452E+00; 8.749611567477869E-02];
vj = [6.071656568115219E-03; -3.994124727031415E-03; -1.194885074562020E-04];
muj = musun/1047.3486;

rs = [-8.602020775180101E+00; -4.532126278101087E+00; 4.209749480415946E-01];
vs = [2.304220049973644E-03; -4.946451041565058E-03; -5.330794761082353E-06];
musa = musun/3497.898;

p0 = [rj,rs,rsun];
v0 = [vj,vs,vsun];
mus = [muj,musa,musun];

%remove any com drift from system
Mtot = sum(mus);
vcom = sum(v0*diag(mus),2)/Mtot;
v0 = v0 - vcom*ones(1,3);

%replace with your favorite integrator
tic;[t,x,dx] = nbodyVect(p0(:),v0(:),mus,0:100*365.25:100000*365.25,'c');toc

r_js = (x(:,1:3) - x(:,7:9)).';
r_ss = (x(:,4:6) - x(:,7:9)).';
v_js = (dx(:,1:3) - dx(:,7:9)).';
v_ss = (dx(:,4:6) - dx(:,7:9)).';

[aj,ej,Ej,Ij,omegaj,Omegaj,Pj,tauj] = vec2orbElem(r_js(:),v_js(:),musun+muj);
[as,es,Es,Is,omegas,Omegas,Ps,taus] = vec2orbElem(r_ss(:),v_ss(:),musun+muj);

%plot
figure(1)
clf
plot3(r_js(1,:),r_js(2,:),r_js(3,:),'.',r_ss(1,:),r_ss(2,:),r_ss(3,:),'.')

%% Secular two-body solution
n1 = sqrt((musun+muj)/aj(1)^3); %rad/day
n2 = sqrt((musun+musa)/as(1)^3);

m = [muj,musa];
n = [n1,n2];
mc = musun;
d = aj(1)/as(1);
b = [3.16269457539456,2.06206350225032];

Apq = @(n,m,d,p,q,b) (-1)^(1 - eq(p,q))*n(p)/4*m(3-p)/(mc+m(p))*d*d^(2-p)*b(2-eq(p,q));

A11 = Apq(n,m,d,1,1,b);
A12 = Apq(n,m,d,1,2,b);
A21 = Apq(n,m,d,2,1,b);
A22 = Apq(n,m,d,2,2,b);

A = [A11, A12; A21, A22];

Bpq = @(n,m,d,p,q,b) (-1)^(1*eq(p,q))*n(p)/4*m(3-p)/(mc+m(p))*d*d^(2-p)*b(1);

B11 = Bpq(n,m,d,1,1,b);
B12 = Bpq(n,m,d,1,2,b);
B21 = Bpq(n,m,d,2,1,b);
B22 = Bpq(n,m,d,2,2,b);

B = [B11, B12; B21, B22];

[fi,lam] = eig(A);
[gi,gam] = eig(B);

%initial conditions
c10 = ej(1)*sin(Omegaj(1)+omegaj(1));
c20 = es(1)*sin(Omegas(1)+omegas(1));
d10 = ej(1)*cos(Omegaj(1)+omegaj(1));
d20 = es(1)*cos(Omegas(1)+omegas(1));
u10 = Ij(1)*sin(Omegaj(1));
u20 = Is(1)*sin(Omegas(1));
v10 = Ij(1)*cos(Omegaj(1));
v20 = Is(1)*cos(Omegas(1));

%eigenvector normalization to match initial conditions
ssa = fi\[c10;c20]; %s1sina1;s2sina2
sca = fi\[d10;d20]; %s1cosa1;s2cosa2
tsb = gi\[u10;u20]; %t1sinb1;t2sinb2
tcb = gi\[v10;v20]; %t1cosb1;t2cosb2

s1s2 = sqrt(ssa.^2 + sca.^2);
a1a2 = atan2(ssa./s1s2,sca./s1s2);

t1t2 = sqrt(tsb.^2 + tcb.^2);
b1b2 = atan2(tsb./t1t2,tcb./t1t2);

%secular time histories
c12 = fi*diag(s1s2)*[sin(lam(1)*t + a1a2(1)).';sin(lam(4)*t + a1a2(2)).'];
d12 = fi*diag(s1s2)*[cos(lam(1)*t + a1a2(1)).';cos(lam(4)*t + a1a2(2)).'];
e12 = sqrt(c12.^2 + d12.^2);

u12 = gi*diag(t1t2)*[sin(gam(1)*t + b1b2(1)).';sin(gam(4)*t + b1b2(2)).'];
v12 = gi*diag(t1t2)*[cos(gam(1)*t + b1b2(1)).';cos(gam(4)*t + b1b2(2)).'];
I12 = sqrt(u12.^2 + v12.^2);

figure(2)
clf
plot(t,ej,'r',t,es,'b',t,e12(1,:),'r--',t,e12(2,:),'b--','Linewidth',2)
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (days)')
ylabel('Eccentricity')
legend({'Jupiter','Saturn'})

figure(3)
clf
plot(t,Ij,'r',t,Is,'b',t,I12(1,:),'r--',t,I12(2,:),'b--','Linewidth',2)
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (days)')
ylabel('Inclination (rad)')
legend({'Jupiter','Saturn'})





