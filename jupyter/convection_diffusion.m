clear all
close all
clc


n=50;
x=linspace(0,10,n);
t=linspace(0,100);
u0=zeros(n,1);
u0(end)=1;

[t,y]=ode45(@PDE_ODEs,t,u0);
pcolor(x,t,y)
shading interp
xlabel('Space')
ylabel('Time')
if max(sum(y,2)-1)<1e-5
title('Accurate to at least 5 d.p.')
end
%%
for i = 1:n
    integral(i) = trapz(x,y(i,:))
end
    

%%
function dydt=PDE_ODEs(t,y)
d=0.01;
c=.5;

middle_terms=d*(y(1:end-2)-2*y(2:end-1)+y(3:end))+c*(y(3:end)-y(2:end-1));

dydt=[d*(y(2)-y(1))+c*y(2);
      middle_terms;
      d*(y(end-1)-y(end))-c*y(end)];

end