function dq=armnew(t,qh)
% q=[1;5;1;0;0;0];
% t=0;
h1=qh(1);
h2=qh(2);
h3=qh(3);
x=qh(4);
y=qh(5);
z=qh(6);
xtx=qh(7);
xty=qh(8);
xtz=qh(9);
h1d=qh(1+9);
h2d=qh(2+9);
h3d=qh(3+9);
xdr=qh(4+9);
ydr=qh(5+9);
zdr=qh(6+9);
xtxdr=qh(7+9);
xtydr=qh(8+9);
xtzdr=qh(9+9);
%取出广义坐标

%% 基本参数设置
h0 = 3;
L0 = 2;
%设置关键参数
dh1=0.00001;
dh2=0.00001;
dh3=0.00001;
dt=0.00001;
h1t=h1+h1d*dt;
h2t=h2+h2d*dt;
h3t=h3+h3d*dt;
m=1.774;
Mp=9.08;
j1=36825.653e-6;
j2=36825.653e-6;
j3=72966.445e-6;
%设置系统参数
M=diag([m,m,m,Mp,Mp,Mp,j1,j2,j3]);
%% 求解变换矩阵部分及其余状态变量部分
Aq=jieA(h1,h2,h3,h0,L0);
Aqh1=jieA(h1+dh1,h2,h3,h0,L0);
Aqh2=jieA(h1,h2+dh2,h3,h0,L0);
Aqh3=jieA(h1,h2,h3+dh3,h0,L0);
At=jieA(h1t,h2t,h3t,h0,L0);
o=Aq*[0;0;0;1];
x=o(1);
y=o(2);
z=o(3);
oh1=Aqh1*[0;0;0;1];
oh2=Aqh2*[0;0;0;1];
oh3=Aqh3*[0;0;0;1];
ot=At*[0;0;0;1];
xd=(ot(1)-x)/dt;
yd=(ot(2)-y)/dt;
zd=(ot(3)-y)/dt;
A=Aq(1:3,1:3);
Az=A';
la=Az*[1;0;0];
lb=Az*[sqrt(3)*L0/2;3*L0/2;h2-h1];
lr=[0;0;1];
la=la./norm(la);
lb=lb./norm(lb);
%求解三个旋转轴的坐标
a2xt=[la,lb,lr];
%角度坐标转换矩阵
%[a,b,r]=jieabr(h1,h2,h3,L0);
[ah1,bh1,rh1]=jieabr(h1+dh1,h2,h3,L0);
[ah2,bh2,rh2]=jieabr(h1,h2+dh2,h3,L0);
[ah3,bh3,rh3]=jieabr(h1,h2,h3+dh3,L0);
[at,bt,rt]=jieabr(h1t,h2t,h3t,L0);
%xt=a2xt*[a;b;r];
xth1=a2xt*[ah1;bh1;rh1];
xth2=a2xt*[ah2;bh2;rh2];
xth3=a2xt*[ah3;bh3;rh3];
xtt=a2xt*[at;bt;rt];
% xtx=xt(1);
% xty=xt(2);
% xtz=xt(3);
xtxd=(xtt(1)-xtx)/dt;
xtyd=(xtt(2)-xty)/dt;
xtzd=(xtt(3)-xtz)/dt;
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
pxtxh1=(xth1(1)-xtx)/dh1;
pxtxh2=(xth2(1)-xtx)/dh2;
pxtxh3=(xth2(1)-xtx)/dh3;
pxtyh1=(xth1(2)-xty)/dh1;
pxtyh2=(xth2(2)-xty)/dh2;
pxtyh3=(xth3(2)-xty)/dh3;
pxtzh1=(xth1(3)-xtz)/dh1;
pxtzh2=(xth2(3)-xtz)/dh2;
pxtzh3=(xth3(3)-xtz)/dh3;
posais=[pxh1,pxh2,pxh3;pyh1,pyh2,pyh3;pzh1,pzh2,pzh3];
posaix=[pxtxh1,pxtxh2,pxtxh3;pxtyh1,pxtyh2,pxtyh3;pxtzh1,pxtzh2,pxtzh3];
posai=[posais,-eye(3),zeros(3,3);posaix,zeros(3,3),-eye(3)];
%% 控制力预设
kg=0;%控制力开关
xm=-0.08;
ym=0.1;
zm=5;
%目标位置
%[h1m,h2m,h3m]=x2hm(xm,ym,zm);
% h1m=2;
% h2m=1;
% h3m=1;
d=posais\[x-xm;y-ym;z-zm];
k=-12.2;%控制力增益
c=18.3;%阻尼力
%F0=kg*[k*(h1-h1m)-c*h1d,k*(h2-h2m)-c*h2d,k*(h3-h3m)-c*h3d];
F0=kg*[k*d(1)-c*h1d,k*d(2)-c*h2d,k*d(3)-c*h3d];
F1=[0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
F0=F0+[0,0.8*sin(pi*t),0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算外力作用
F1s=Az*[0;0;-F1(1)-F1(2)-F1(3)];
F1j=[F1(3)*L0;sqrt(3)*F1(2)*L0/2-sqrt(3)*F1(1)*L0/2;0];
Q=[F0(1);F0(2);F0(3);F1s;F1j];
%% 动力学矩阵正式得出
q=[h1d;h2d;h3d;xdr;ydr;zdr;xtxdr;xtydr;xtzdr];
%kx=jiekx(h1t,h2t,h3t,dh1,dh2,dh3,dt,posai,h0,L0,q);
kx=jiekx0new(h1,h2,h3,h1t,h2t,h3t,dh1,dh2,dh3,dt,q);
qf=[Q;kx];
ABl=[M,posai';posai,zeros(6,6)];
niABl=inv(ABl)/2;
%ABr=[zeros(9,6),eye(9),zeros(9,9);zeros(9,15),eye(9);zeros(6,24)];
dql=niABl*qf;%带拉格朗日乘子的状态矩阵
dq=[h1d;h2d;h3d;xdr;ydr;zdr;xtxdr;xtydr;xtzdr;dql(1);dql(2);dql(3);dql(4);dql(5);dql(6);dql(7);dql(8);dql(9)];%只含h的状态矩阵
end