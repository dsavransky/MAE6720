%% An animation of the eliptic flyaround problem 
%%
load('topo.mat','topo','topomap1');

%%
R_E = 6378.137; %km
h = 420; %km
a = h + R_E; 
mu = 398600.4418; %km^2/s^3
n = sqrt(mu/a^3);
y0 = -25;

%flyaround
tf = 2*pi/n;
dx0 = -y0*n/2;

t = linspace(0,tf,100);
xss = dx0*sin(n*t)/n;
yss = 2*dx0*(cos(n*t) - 1)/n;
xss = cos(n*t).*xss - sin(n*t).*yss;
yss = sin(n*t).*xss + cos(n*t).*yss;
%%
figure(1)
clf

[x,y,z] = sphere(50);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
Earth = surface(x,y,z,props);

d = 1.5;

props2.AmbientStrength = 0.1;
props2.DiffuseStrength = 1;
props2.SpecularColorReflectance = .5;
props2.SpecularExponent = 20;
props2.SpecularStrength = 1;
props2.FaceColor= 'blue';
props2.EdgeColor = 'none';
iss = surface(x/20+d,y/20,z/20,props2);


props2.FaceColor = 'red';
ss = surface(x/20+d,y/20-0.2,z/20,props2);
%rotate(ss,[1,0,0],60,[d*cos(th0)+0.1,d*sin(th0)+0.1,0])

% Add lights.
l = light('position',[2 0 0.5]);
view(35,10)
axis equal off

%%
xss1 = xss*0.2/25;
yss1 = yss*0.2/25;
for j = 2:length(t)
    set(ss,'XData',get(ss,'XData')+xss1(j)-xss1(j-1),'YData',get(ss,'YData')+yss1(j)-yss1(j-1))
      
    xs = get(Earth,'XData');
    ys = get(Earth,'YData');
    zs = get(Earth,'ZData');
    vns = get(Earth,'VertexNormals');
    n = length(xs);
    
    th = 2*pi/length(t);
    rotmat = [cos(th), -sin(th),0;sin(th), cos(th),0;0,0,1];
    
    tmp = rotmat*[xs(:),ys(:),zs(:)].';
    newx = reshape(tmp(1,:),n,n);
    newy = reshape(tmp(2,:),n,n);
    newz = reshape(tmp(3,:),n,n);
    
    tmp2 = reshape(vns,n*n,3).';
    newvns = rotmat*tmp2;
    newvns = reshape(newvns.',n,n,3);
    
    set(Earth,'XData',newx,'YData',newy,'ZData',newz,'VertexNormals',newvns)
    set(l,'position',[2*cos(2*pi*t(j)/tf),2*sin(2*pi*t(j)/tf),0.5])
    pause(0.1)
end

%%
figure(2)
clf

[x,y,z] = sphere(50);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
Earth2 = surface(x,y,z,props);

d = 1.5;

props2.AmbientStrength = 0.1;
props2.DiffuseStrength = 1;
props2.SpecularColorReflectance = .5;
props2.SpecularExponent = 20;
props2.SpecularStrength = 1;
props2.FaceColor= 'blue';
props2.EdgeColor = 'none';
iss2 = surface(x/20+d,y/20,z/20,props2);

axis equal off
%%
xss2 = xss*0.15/25;
yss2 = yss*0.15/25;

set(gca,'CameraPosition',[d,-0.15,0])
set(gca,'CameraTarget',[d,0,0])
set(gca,'CameraViewAngle',140)


for j = 2:length(t)
    set(gca,'CameraPosition',[d+xss2(j),-0.2+yss2(j),0])

    xs = get(Earth,'XData');
    ys = get(Earth,'YData');
    zs = get(Earth,'ZData');
    vns = get(Earth,'VertexNormals');
    n = length(xs);
    
    th = 2*pi/length(t);
    rotmat = [cos(th), -sin(th),0;sin(th), cos(th),0;0,0,1];
    
    tmp = rotmat*[xs(:),ys(:),zs(:)].';
    newx = reshape(tmp(1,:),n,n);
    newy = reshape(tmp(2,:),n,n);
    newz = reshape(tmp(3,:),n,n);
    
    tmp2 = reshape(vns,n*n,3).';
    newvns = rotmat*tmp2;
    newvns = reshape(newvns.',n,n,3);
    
    set(Earth2,'XData',newx,'YData',newy,'ZData',newz,'VertexNormals',newvns)
    pause(0.1)
end