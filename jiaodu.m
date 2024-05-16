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

t=0:0.1:15;
hh1=1+sin(t-120*pi/180);
hh2=1+sin(t+0*pi/180);
hh3=1+sin(t+120*pi/180);
xit=[];
L0 = 2;
h0 = 3;

for i=1:length(hh3)

    h1 = 3+hh1(i);
    h2 = 3+hh2(i);
    h3 = 3+hh3(i);
    xita2=[(3*(h3-h1)*L0)/2,1.732*((h2-h1)*L0-((h3-h1)*L0)/2),-3*1.732*L0*L0/2];
    aa=xita2;
    chui=[0,0,1];
    xx1=acos(dot( aa,chui)/(norm( aa)*norm(chui)));
    L13f = 3*L0*L0+(h3-h1)*(h3-h1);
L23f = 3*L0*L0+(h3-h2)*(h3-h2);  %L13f即L13的平方

    L1 = sqrt(9*L0*L0/4+(h2-h1)*(h2-h1))-L0/2;
    [r,xx,yy] = gammaqiujie(L0,L1,L13f,L23f);
    b = cosbeta1hanshu(h1,h2,h3,L0);
    alpha = atan((h2-h1)/(3*L0/2));
    Lq = (1.732*L0*L1/2)/sqrt(L0*L0+L1*L1+L0*L1);
    xit=[xit;b];

end
plot(t,xit.*180./pi);
%
% plot(t,x)
% hold on
% plot(t,y)
% plot(t,z)
% hold off


