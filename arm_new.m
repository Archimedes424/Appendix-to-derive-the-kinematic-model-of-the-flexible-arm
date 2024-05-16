function dq=arm_new(t,qh)
% q=[1;5;1;0;0;0];
% t=0;
h1=qh(1);
h2=qh(2);
h3=qh(3);
h1d=qh(4);
h2d=qh(5);
h3d=qh(6);
%取出广义坐标

%% 基本参数设置
h0 = 350e-3;
L0 = 55e-3;
%设置关键参数
dh1=0.000001;
dh2=0.000001;
dh3=0.000001;
dt=0.000001;
h1t=h1+h1d*dt;
h2t=h2+h2d*dt;
h3t=h3+h3d*dt;
m=1.774;
Mp=9.08;
j1=36825.653e-6;
j2=36825.653e-6;
j3=72966.445e-6;
%设置系统参数
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
[a,b,r]=jieabr(h1,h2,h3,L0);
[ah1,bh1,rh1]=jieabr(h1+dh1,h2,h3,L0);
[ah2,bh2,rh2]=jieabr(h1,h2+dh2,h3,L0);
[ah3,bh3,rh3]=jieabr(h1,h2,h3+dh3,L0);
[at,bt,rt]=jieabr(h1t,h2t,h3t,L0);
ad=(at-a)/dt;
bd=(bt-b)/dt;
rd=(rt-r)/dt;
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
F0=F0+[0.1045*sin(pi*t),0,0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算外力作用
F1s=Az*[0;0;-F1(1)-F1(2)-F1(3)];
F1j=[F1(3)*L0;sqrt(3)*F1(2)*L0/2-sqrt(3)*F1(1)*L0/2;0];
Q=[F0(1);F0(2);F0(3);F1s;F1j];
%% 动力学矩阵正式得出
M=[diag([m,m,m,Mp,Mp,Mp]),zeros(6,3);zeros(3,6),a2xt'*diag([j1,j2,j3])*a2xt];
q=[h1d;h2d;h3d;xd;yd;zd;ad;bd;rd];
%kx=jiekx(h1t,h2t,h3t,dh1,dh2,dh3,dt,posai,h0,L0,q);
kx=jiekx0new11(h1,h2,h3,h1t,h2t,h3t,dh1,dh2,dh3,dt,q);
qf=[Q;kx];
ABl=[M,posai';posai,zeros(6,6)];
%ABr=[zeros(9,6),eye(9),zeros(9,9);zeros(9,15),eye(9);zeros(6,24)];
dql=ABl\qf;%带拉格朗日乘子的状态矩阵
dq=[h1d;h2d;h3d;dql(1);dql(2);dql(3)];%只含h的状态矩阵
end