close all
clear all

n = 100;
x = linspace(0,3000,n);
t = linspace(0,3600*8,n);
m = 1;
D = 100;
sol = pdepe(m,@diffusion,@diffusionic,@diffusionbc,x,t);
u = sol(:,:,1);

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


%----------------------------------------------
function [c,f,s] = diffusion(x,t,u,dudx)
D = 100;
c = 1;
r = - 0.00001;
f = D*dudx; %fixed domain
s = 0; %fixed domain
end
%----------------------------------------------
function u0 = diffusionic(x)
D = 100;
T = 1;
u0 = 1/(4*pi*D*T)*exp(-x^2/(4*D*T));
end
%----------------------------------------------
function [pl,ql,pr,qr] = diffusionbc(xl,ul,xr,ur,t)
pl = 0; %ignored by solver since m=1
ql = 1; %ignored by solver since m=1
pr = 0; 
qr = 1;
end
%-----------------------------------------------
function one = Intone(x,t,u,n)
    one = t;
    for i = 1:n
        one(i) = trapz(x,u(i,:).*x.*2.*pi);
    end
end
%------------------------------------------------
function r2 = Intr2(x,t,u,n)
    r2 = t;
    for i = 1:n
        r2(i) = trapz(x,u(i,:).*x.^2.*2.*pi);
    end
end