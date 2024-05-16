function posaiqqd=pqqd11(h1,h2,h3,h1d,h2d,h3d,dh1,dh2,dh3,dt,qd)
h0 = 3;
L0 = 2;
Aq=jieA(h1,h2,h3,h0,L0);
Aqh1=jieA(h1+dh1,h2,h3,h0,L0);
Aqh2=jieA(h1,h2+dh2,h3,h0,L0);
Aqh3=jieA(h1,h2,h3+dh3,h0,L0);
o=Aq*[0;0;0;1];
x=o(1);
y=o(2);
z=o(3);
oh1=Aqh1*[0;0;0;1];
oh2=Aqh2*[0;0;0;1];
oh3=Aqh3*[0;0;0;1];
A=Aq(1:3,1:3);
Az=A';
la=Az*[1;0;0];
lb=Az*[sqrt(3)*L0/2;3*L0/2;h2-h1];
lr=[0;0;1];
la=la./norm(la);
lb=lb/norm(lb);
%求解三个旋转轴的坐标
a2xt=[la,lb,lr];
%角度坐标转换矩阵
[a,b,r]=jieabr(h1,h2,h3,L0);
[ah1,bh1,rh1]=jieabr(h1+dh1,h2,h3,L0);
[ah2,bh2,rh2]=jieabr(h1,h2+dh2,h3,L0);
[ah3,bh3,rh3]=jieabr(h1,h2,h3+dh3,L0);

%% 求数值偏微分部分
pxh1=(oh1(1)-x)/dh1;
pxh2=(oh2(1)-x)/dh2;
pxh3=(oh3(1)-x)/dh3;
pyh1=(oh1(2)-y)/dh1;
pyh2=(oh2(2)-y)/dh2;
pyh3=(oh3(2)-y)/dh3;
pzh1=(oh1(3)-z)/dh1;
pzh2=(oh2(3)-z)/dh2;
pzh3=(oh3(3)-z)/dh3;
pah1=(ah1-a)/dh1;
pah2=(ah2-a)/dh2;
pah3=(ah3-a)/dh3;
pbh1=(bh1-b)/dh1;
pbh2=(bh2-b)/dh2;
pbh3=(bh3-b)/dh3;
prh1=(rh1-r)/dh1;
prh2=(rh2-r)/dh2;
prh3=(rh3-r)/dh3;
posais=[pxh1,pxh2,pxh3;pyh1,pyh2,pyh3;pzh1,pzh2,pzh3];
posaix=[pah1,pah2,pah3;pbh1,pbh2,pbh3;prh1,prh2,prh3];
posai=[posais,-eye(3),zeros(3,3);posaix,zeros(3,3),-eye(3)];
posaiqqd=posai*qd;
end