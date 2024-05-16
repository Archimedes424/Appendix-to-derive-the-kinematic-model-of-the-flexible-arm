close all
clear
clc
%读取仿真数据
data=readmatrix("动力学\新建 XLSX 工作表(2).xlsx");
tf=data(:,1);
yf=data(:,2);
xf=data(:,3);
zf=data(:,4);
h1f=data(:,5);
h2f=data(:,6);
h3f=data(:,7);
q0=[350e-3;350e-3;350e-3;0;0;0];
%tspan=[0;20];
tspan=tf;
[t,h]=ode45(@arm_new,tspan,q0);
k=1e3;
h1=k*h(:,1);
h2=k*h(:,2);
h3=k*h(:,3);
[x,y,z,xtx,xty,xtz,xd,yd,zd,xtxd,xtyd,xtzd]=quan(h);
x=k*x;
y=k*y;
z=k*z;
subplot(3,2,1)
plot(t,h1,Color='b',LineWidth=2);
hold on
plot(tf,h1f,LineStyle='none',Marker='o',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('h1,1(mm)')
legend('Kinematic equations','Software emulation')
subplot(3,2,3)
plot(t,h2,Color='b',LineWidth=2,LineStyle='--');
hold on
plot(tf,h2f,LineStyle='none',Marker='o',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('h1,2(mm)')
legend('Kinematic equations','Software emulation')
subplot(3,2,5)
plot(t,h3,Color='b',LineWidth=2,LineStyle='-.');
hold on
plot(tf,h3f,LineStyle='none',Marker='o',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('h1,3(mm)')
legend('Kinematic equations','Software emulation')

subplot(3,2,2)
plot(t,x,Color='r',LineWidth=1.7);
hold on
plot(tf,xf,LineStyle='none',Marker='^',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('x(mm)')
legend('Kinematic equations','Software emulation')
subplot(3,2,4)
plot(t,y,Color='r',LineWidth=1.7,LineStyle='--');
hold on
plot(tf,yf,LineStyle='none',Marker='^',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('y(mm)')
legend('Kinematic equations','Software emulation')
subplot(3,2,6)
plot(t,z,Color='r',LineWidth=1.7,LineStyle='-.');
hold on
plot(tf,zf,LineStyle='none',Marker='^',MarkerSize=5,Color='k')
hold off
xlabel('t(s)')
ylabel('z(mm)')
legend('Kinematic equations','Software emulation')