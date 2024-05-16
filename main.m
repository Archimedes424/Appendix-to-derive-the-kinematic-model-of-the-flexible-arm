close all
clear
clc
%% 求解部分
q0=[350e-3;350e-3;350e-3;0;0;0];
%tspan=[0;20];
tspan=0:0.01:10;
[t,h]=ode45(@arm_new,tspan,q0);
%% 绘图数据计算部分
%反算出所有坐标的变换方输
h1=h(:,1);
h2=h(:,2);
h3=h(:,3);
h1d=h(:,4);
h2d=h(:,5);
h3d=h(:,6);
[x,y,z,xtx,xty,xtz,xd,yd,zd,xtxd,xtyd,xtzd]=quan(h);
%% 三根杆位置的变化图
subplot(3,2,1)
plot(t,h1*1e3,'LineWidth',2)
title('第1根杆的位置变化(mm)')
subplot(3,2,3)
plot(t,h2*1e3,'LineWidth',2)
title('第2根杆的位置变化(mm)')
subplot(3,2,5)
plot(t,h3*1e3,'LineWidth',2)
title('第3根杆的位置变化(mm)')
%% 盘心三个坐标变化图
subplot(3,2,2)
plot(t,x*1e3,'r','LineWidth',2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
title('x的位置变化(mm)')
subplot(3,2,4)
plot(t,y*1e3,'r','LineWidth',2)
title('y的位置变化(mm)')
subplot(3,2,6)
plot(t,z*1e3,'r','LineWidth',2)
title('z的位置变化(mm)')
%% 盘心位置的变化图
% scatter3(x,y,z,'b')
% title('盘心位置的变化图')