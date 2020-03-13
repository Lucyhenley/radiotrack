close all
clear all

r =  - 0.005;
n = 100;
x = linspace(0,2000,n);
t = linspace(0,25000,n);
m = 0;
D = 75;
T = 4000;

sol = pdepe(m,@diffusion,@diffusionic,@diffusionbc,x,t);

u = sol(:,:,1);
T=4000;

figure
s=surf(x,t,u);
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
colorbar
s.EdgeColor='None';
view(2)

figure
s = plot(x,u(1,:)); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 
v = VideoWriter('u_apical.avi');
v.FrameRate = 5
open(v);
for k = 1:n 
   plot(x,u(k,:))
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);

figure
one = moment(x,t,u,n,m,0);
scatter(t,one)
xlabel('t')
ylabel('$\int_0^R \phi dx$','Interpreter', 'latex')

figure
r2 = moment(x,t,u,n,m,2);
scatter(t,r2)
xlabel('t')
ylabel('MSD')

function u0 = IC(x)
    D = 75;
    T = 4000;
    u0 = 1/(4*pi*D*T)*exp(-x^2/(4*D*T));
end
%---------------------------------
function one = moment(x,t,u,n,m,O)
    one = t;
    for i = 1:n
        one(i) = trapz(x,u(i,:).*x.^(m+O).*(2*pi)^m);
    end
end
%----------------------------------------------
function l = lt(r,t)
    l = 1 + r*t^2; 
    %l = exp(r*t);
end
%---------------------------------
function dl = dlt(r,t)
    dl = 2*r*t;
   % dl = r*exp(r*t);
end
%---------------------------------
function [c,f,s] = diffusion(x,t,u,dudx)
D = 10;
c = 1;
r =  - 0.00005;
l = lt(r,t);
dl = dlt(r,t);
%f = D*dudx; %fixed domain
f = D*dudx/l^2 + dl*x*u/l;% %apical growth
%f = D*dudx/l^2; % isotropic growth
%s = 0; %fixed domain
s = -r*u ; %apical growth
%s = - dl*u/l; %isotropic growth
end
%----------------------------------------------
function u0 = diffusionic(x)
D = 75;
T = 4000;
%u0 = 1/(4*pi*D*T)*exp(-x^2/(4*D*T));
u0 = 2/sqrt(2*pi*D*T)*exp(-x^2/(2*D*T));
end
%----------------------------------------------
function [pl,ql,pr,qr] = diffusionbc(xl,ul,xr,ur,t)
r = -0.001;
pl = 0; 
ql = 1; 
pr = 0;
qr = 1; 
end
%---------------------