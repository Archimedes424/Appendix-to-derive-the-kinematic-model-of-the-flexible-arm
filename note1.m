close all
clear
clc
% hh1 = 'h1?  ';
% hh2 = 'h2?  ';
% hh3 = 'h3?  ';
% hh0 = 'h0?  ';
% LL0 = 'L0?  ';
% 
% h0 = input(hh0);
% h1 = input(hh1);
% h2 = input(hh2);
% h3 = input(hh3);
% L0 = input(LL0);

t=0:0.01:10;
hh1=50*sin(pi*t/2+30*pi/180);
hh2=50*sin(2*pi*t/2+70*pi/180);
hh3=50*sin(pi*t/2+90*pi/180);
x=[];
y=[];
z=[];
for i=1:length(hh3)
h0 = 349.5;
h1 = h0+hh1(i);
h2 = h0+hh2(i);
h3 = h0+hh3(i);
L0 = 60;

A1 = [1 0 0 0;
      0 1 0 0;
      0 0 1 h0;
      0 0 0 1];
A2 = [1 0 0 0;
      0 1 0 0;
      0 0 1 h1-h0;
      0 0 0 1];
alpha = atan((h2-h1)/(3*L0/2));
A3 = [1,0,0,0;
      0,cos(alpha),-sin(alpha),-L0/2+L0*cos(alpha)/2;
      0,sin(alpha),cos(alpha),L0*sin(alpha)/2;
      0,0,0,1;];

%L1 = sqrt(3*L0*L0+(h2-h1)*(h2-h1))-L0/2;这个好像不太对 应该是3/2？
L1 = sqrt(9*L0*L0/4+(h2-h1)*(h2-h1))-L0/2;

b = cosbeta1hanshu(h1,h2,h3,L0);

Lq = (sqrt(3)*L0*L1/2)/sqrt(L0*L0+L1*L1+L0*L1);

x1 = -Lq*(1-cos(b))*((L1+L0/2)/sqrt(L0*L0+L1*L1+L0*L1));
x2 = Lq*(1-cos(b))*((sqrt(3)*L0/2)/sqrt(L0*L0+L1*L1+L0*L1));
x3 = -Lq*sin(b);
zhu=[sqrt(3)*L0/2,L0/2+L1,0];
p1=zhu(1)/sqrt(sum(zhu.^2));
p2=zhu(2)/sqrt(sum(zhu.^2));
A4 = [p1*p1*(1-cos(b))+cos(b),p1*p2*(1-cos(b)),p2*sin(b),x1;
      p1*p2*(1-cos(b)),p2*p2*(1-cos(b))+cos(b),-p1*sin(b),x2;
      -p2*sin(b),p1*sin(b),cos(b),x3;
      0,0,0,1;];

r0 = sqrt(L0*L0+L0*L1+L1*L1)/sqrt(3);  % r0即为新的圆半径
L13f = 3*L0*L0+(h3-h1)*(h3-h1);
L23f = 3*L0*L0+(h3-h2)*(h3-h2);  %L13f即L13的平方
 
%[r,xx,yy] = gammaqiujie(L0,L1,L13f,L23f);
[r,xx,yy] = gammafast(h1,h2,h3,L13f,L23f,L0,L1);
% r=0;
% xx=0;
% yy=0;

A5=[cos(r) -sin(r) 0 xx;
    sin(r) cos(r) 0 yy;
    0 0 1 0;
    0 0 0 1];

T = A1*A2*A3*A4*A5;
x3=T*[0;0;0;1];
x=[x;x3(1)];
y=[y;x3(2)];
z=[z;x3(3)];
end
subplot(3,1,1)
plot(t,x)
title('x')
subplot(3,1,2)
plot(t,y);
title('y')
subplot(3,1,3)
plot(t,z);
title('z')
% 
% plot(t,x)
% hold on
% plot(t,y)
% plot(t,z)
% hold off


