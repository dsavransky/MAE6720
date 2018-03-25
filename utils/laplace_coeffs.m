% A variety of ways to calculate the laplace coeficient and its derivative

%via the hypergeometric function
lc = @(s,k,alpha) 2*pochhammer(s,abs(k))*hypergeom([s,s+abs(k)],abs(k)+1,alpha^2)*alpha^abs(k)/factorial(abs(k));

%via direct quadrature 
fun = @(x,s,k,alpha) cos(k*x)./(1 - 2*alpha*cos(x) + alpha.^2).^s;
lc2 = @(s,k,alpha) integral(@(x)fun(x,s,k,alpha),0,2*pi)/pi;

%alternative quadrature formulation
fun2 = @(x,s,k,h,alpha) cos(k*x).*(cos(x) - alpha).^h./(1 - 2*alpha*cos(x) + alpha.^2).^s;
nhs =  @(s,k,h,alpha) integral(@(x)fun2(x,s,k,h,alpha),0,pi)*2/pi;
lc3 = @(s,k,alpha) nhs(s,k,0,alpha);

%first and second derivatives via recursion equations
Dlc = @(s,k,alpha) s*(lc(s+1,k-1,alpha) - 2*alpha*lc(s+1,k,alpha) + lc(s+1,k+1,alpha));
D2lc = @(s,k,alpha) s*(Dlc(s+1,k-1,alpha) - 2*alpha*Dlc(s+1,k,alpha) + Dlc(s+1,k+1,alpha) - 2*lc(s+1,k,alpha))

%first and second derivatives by quadrature
Dlc2 = @(s,k,alpha) 2*s*nhs(s+1,k,1,alpha);
D2lc2 = @(s,k,alpha) 4*s*(s+1)*nhs(s+2,k,2,alpha) - 2*s*nhs(s+1,k,0,alpha);