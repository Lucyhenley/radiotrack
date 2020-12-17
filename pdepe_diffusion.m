close all
clear all

n = 1000;
nt = 10000;
R = 2000;
x = linspace(0,R,n);
t = linspace(0,3600*8,nt);
m = 1;
D = 100;
sol = pdepe(m,@diffusion,@diffusionic,@diffusionbc,x,t);
u = sol(:,:,1);

figure
one = Intone(x,t,u,nt);
plot(t,one)
xlabel('t')
ylabel('Integral ux dx * 2pi')

figure
plot(t,u(:,end))
xlabel('t')
ylabel('u(R)')

figure
r2 = Intr2(x,t,u,nt);
r2calc = r2calcs(x,t,u,nt,D,R);
plot(t,r2)
hold on
plot(t,r2calc)
yline(R^2/2)
legend({'Integral r^2 \phi rdr d\theta','Integral solution','R^2/2'},'Location','southeast')
xlabel('t')
ylabel('Integral')


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
T = 100;
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
function one = Intone(x,t,u,n,dx)
    one = t; 
    for i = 1:n
        one(i) = trapz(x,u(i,:).*x.*2.*pi);
    end
end
%------------------------------------------------
function r2 = Intr2(x,t,u,n)
    r2 = t;
    for i = 1:n
        r2(i) = trapz(x,u(i,:).*x.^3.*2.*pi);
    end
end
%------------------------------------------------
%4D*(ts - ?*L^2 * trapz(uR2d[1:i],t[1:i]))
function r2calc = r2calcs(x,t,u,nt,D,R)
    r2calc = t;
    for i = 2:length(t)
        r2calc(i) = 4*D*(t(i) - pi*R^2*(trapz(t(1:i),u(1:i,end))));
    end
end