close all
clear all

n = 100;
x = linspace(0,3000,n);
t = linspace(0,3600*8,n);
m = 1;
D = 75;
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
one = Intone(x,t,u,n);
scatter(t,one)
xlabel('t')
ylabel('Integral ux dx * 2pi')

figure
r2 = Intr2(x,t,u,n);
scatter(t,r2)
xlabel('t')
ylabel('Integral ux dx * 2pi * x^2')

function u0 = IC(x)
    D = 75;
    T = 4000;
    u0 = 1/(4*pi*D*T)*exp(-x^2/(4*D*T));
    
end
%---------------------------------
function one = Intone(x,t,u,n)
    one = t;
    for i = 1:n
        one(i) = trapz(x,u(i,:).*x.*2.*pi);
    end
end
%---------------------------------
function r2 = Intr2(x,t,u,n)
    r2 = t;
    for i = 1:n
        r2(i) = trapz(x,u(i,:).*x.^2.*2.*pi);
    end
end
%----------------------------------------------
function l = lt(r,t)
    l = 1500*exp(r*t);
end
%---------------------------------
function dl = dlt(r,t)
    dl = 1500*r*exp(r*t);
end
%---------------------------------
function [c,f,s] = diffusion(x,t,u,dudx)
l = -0.0001;
D = 75;
c = 1;
r = - 0.00001;
f = D*dudx; %fixed domain
f = D*dudx/exp(2*r*t);  %dlt(r,t)/lt(r,t); %growing domain
s = 0; %fixed domain
s = r*dudx; %growing domain 
end
%----------------------------------------------
function u0 = diffusionic(x)
D = 75;
T = 4000;
u0 = 1/(4*pi*D*T)*exp(-x^2/(4*D*T));
end
%----------------------------------------------
function [pl,ql,pr,qr] = diffusionbc(xl,ul,xr,ur,t)
r = -0.00001;
dr = 1;
D = 75;
pl = 0; %ignored by solver since m=1
ql = 1; %ignored by solver since m=1
pr = 0; 
qr = 1;
pr = -r*xr*ur;
qr = exp(2*r*t)/D;
end
%---------------------