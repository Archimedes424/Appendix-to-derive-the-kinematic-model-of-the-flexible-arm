function ykb=suanykb(h1,h2,h3,h0,L0)
dh1=0.0001;
dh2=0.0001;
dh3=0.0001;
%% 求解变换矩阵部分及其余状态变量部分
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
ykb=[pxh1,pxh2,pxh3;pyh1,pyh2,pyh3;pzh1,pzh2,pzh3];
end