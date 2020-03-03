l = 1:4;
nu = 1:4;

t = linspace(0,2*pi,100);
figure(1)
clf
counter = 1;
for j = 1:length(l)
    for k = 1:length(nu)
        subplot(length(l),length(nu),counter)
        x = cos(l(j)*t);
        y = sin(l(j)*t);
        z = sin(nu(k)*t);
        plot3(x,y,z)
        set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'FontSize',16)
        if j == 1
            title(['\nu = ',num2str(nu(k))])
        end
        if k == 1
            zlabel(['\lambda = ',num2str(l(j))])
        end
        counter = counter + 1;
        
    end
end