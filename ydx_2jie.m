close all
clear
clc
t=0:0.01:5;%运动时间范围
kh1v1=50*sin(pi*t+0*pi/180);
kh1v2=50*sin(2*pi*t+120*pi/180);
kh1v3=50*sin(pi*t-120*pi/180);%下盘面三杆伸缩情况
kh2v1=50*sin(2*pi*t-120*pi/180);
kh2v2=50*sin(2*pi*t+0*pi/180);
kh2v3=50*sin(2*pi*t-120*pi/180);%上盘面三杆伸缩情况
%长度单位均为mm
L0 = 60;
h0=350;
x1=[];
y1=[];
z1=[];
x2=[];
y2=[];
z2=[];
z=[-1,0,0,0;0,-1,0,0;0,0,1,0;0,0,0,1];%中转矩阵
%% 运动学解算
for i=1:length(t)
h1v1=h0+kh1v1(i);
h1v2=h0+kh1v2(i);
h1v3=h0+kh1v3(i);
h2v1=h0+kh2v1(i)+5;
h2v2=h0+kh2v2(i)+5;
h2v3=h0+kh2v3(i)+5;
A1=jieA(h1v1,h1v2,h1v3,h0,L0);
A2=jieA(h2v1,h2v2,h2v3,h0,L0);
r1=A1*[0;0;0;1];
r2=A1*z*A2*[0;0;0;1];
x1=[x1;r1(1)];
y1=[y1;r1(2)];
z1=[z1;r1(3)];
x2=[x2;r2(1)];
y2=[y2;r2(2)];
z2=[z2;r2(3)];
end
%% 后处理
subplot(1,2,1)
title('盘面1质心运动')
subplot(1,2,2)
title('盘面2质心运动')
subplot(3,2,1)
plot(t,x1,LineWidth=1,Color='r');
xlabel('t(s)')
ylabel('x1')
subplot(3,2,3)
plot(t,y1,LineWidth=1,Color='b');
xlabel('t(s)')
ylabel('y1')
subplot(3,2,5)
plot(t,z1,LineWidth=1,Color='m')
xlabel('t(s)')
ylabel('z1')
subplot(3,2,2)
plot(t,x2,LineWidth=1,Color='r',LineStyle=':');
xlabel('t(s)')
ylabel('x2')
subplot(3,2,4)
plot(t,y2,LineWidth=1,Color='b',LineStyle=':');
xlabel('t')
ylabel('y2')
subplot(3,2,6)
plot(t,z2,LineWidth=1,Color='m',LineStyle=':');
xlabel('t(s)')
ylabel('z2')
%三维轨迹图
figure
scatter3(x1,y1,z1)
hold on
scatter3(x2,y2,z2)
%% 数据对比
data=readmatrix("数据.xlsx");
tc=data(:,1);
x2c=data(:,2);
y2c=data(:,3);
z2c=data(:,4);
x1c=data(:,5);
y1c=data(:,6);
z1c=data(:,7);

subplot(1,2,1)
title('盘面1质心运动')
subplot(1,2,2)
title('盘面2质心运动')
subplot(3,2,1)
plot(t,x1,LineWidth=2,Color='r');
xlabel('t(s)')
ylabel('x1(mm)')
hold on
plot(tc,x1c,LineWidth=0.6,Color='k',LineStyle='none',Marker='^');
hold off
legend('Kinematic equations','Software emulation')
subplot(3,2,3)
plot(t,y1,LineWidth=2,Color='b');
xlabel('t(s)')
ylabel('y1(mm)')
hold on
plot(tc,y1c,LineWidth=0.6,Color='k',LineStyle='none',Marker='*');
hold off
legend('Kinematic equations','Software emulation')
subplot(3,2,5)
plot(t,z1,LineWidth=2,Color='m')
xlabel('t(s)')
ylabel('z1(mm)')
hold on
plot(tc,z1c,LineWidth=0.6,Color='k',LineStyle='none',Marker='+');
hold off
legend('Kinematic equations','Software emulation')
subplot(3,2,2)
plot(t,x2,LineWidth=2,Color='r',LineStyle=':');
xlabel('t(s)')
ylabel('x2(mm)')
hold on
plot(tc,x2c,LineWidth=0.6,Color='k',LineStyle='none',Marker='^');
hold off
legend('Kinematic equations','Software emulation')
subplot(3,2,4)
plot(t,y2,LineWidth=2,Color='b',LineStyle=':');
xlabel('t(s)')
ylabel('y2(mm)')
hold on
plot(tc,y2c,LineWidth=0.6,Color='k',LineStyle='none',Marker='*');
hold off
legend('Kinematic equations','Software emulation')
subplot(3,2,6)
plot(t,z2,LineWidth=2,Color='m',LineStyle=':');
xlabel('t(s)')
ylabel('z2(mm)')
hold on
plot(tc,z2c,LineWidth=0.6,Color='k',LineStyle='none',Marker='+');
hold off
legend('Kinematic equations','Software emulation')
