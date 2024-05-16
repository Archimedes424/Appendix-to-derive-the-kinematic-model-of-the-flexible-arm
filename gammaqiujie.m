function [r,x,y] = gammaqiujie(L0,L1,L13f,L23f)
%P1P3=(x13,y13)
syms x13 y13;
f3 = x13*x13+y13*y13-L13f;
f4 = (x13-sqrt(3)*L0/2)*(x13-sqrt(3)*L0/2)+(L0/2+L1-y13)*(L0/2+L1-y13)-L23f;
p13mo=L13f;
p23mo=L23f;
l0=L0;
l1=L1;
[x13,y13] = solve(f3,f4,x13>0);
x13 = double(x13);
y13 = double(y13); 
% x1301=-(-3*l0^4 - 3*l0^3*l1 - 3*l0^2*l1^2 - 3*l0^2*p13mo^2 + 3*l0^2*p23mo^2 + sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2))*l0 + 2*sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2))*l1)*sqrt(3)/(12*(l0^2 + l0*l1 + l1^2)*l0);
% x1302=(3*l0^4 + 3*l0^3*l1 + 3*l0^2*l1^2 + 3*l0^2*p13mo^2 - 3*l0^2*p23mo^2 + sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2))*l0 + 2*sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2))*l1)*sqrt(3)/(12*(l0^2 + l0*l1 + l1^2)*l0);
% y1301=-(-l0^3 - 3*l0^2*l1 - 3*l0*l1^2 - l0*p13mo^2 + l0*p23mo^2 - 2*l1^3 - 2*l1*p13mo^2 + 2*l1*p23mo^2 + sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2)))/(4*(l0^2 + l0*l1 + l1^2));
% y1302=(l0^3 + 3*l0^2*l1 + 3*l0*l1^2 + l0*p13mo^2 - l0*p23mo^2 + 2*l1^3 + 2*l1*p13mo^2 - 2*l1*p23mo^2 + sqrt(-3*l0^2*(l0^2 + l0*l1 + l1^2 - p13mo^2 - 2*p13mo*p23mo - p23mo^2)*(l0^2 + l0*l1 + l1^2 - p13mo^2 + 2*p13mo*p23mo - p23mo^2)))/(4*(l0^2 + l0*l1 + l1^2));
% x13=max(x1301,x1302);
% y13=min(y1301,y1302);

%//在O3系中用120度关系求解Oc'！r0好像没啥用啊
syms x y;
OcpP1 = [-sqrt(3)*L0/2-x,-L0/2-y];
P1P3 = [x13,y13];
P1P2 = [sqrt(3)*L0/2,L0/2+L1]; % p.s.全用的横着的！
f1 = (sqrt(3)*L0/3+sqrt(3)*L1/6)*(sqrt(3)*L0/3+sqrt(3)*L1/6)+(-L1/2)*(-L1/2)-(sqrt(3)*L0/3+sqrt(3)*L1/6+x)*(sqrt(3)*L0/3+sqrt(3)*L1/6+x)-(-L1/2+y)*(-L1/2+y);
f2 = norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
f3=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);

[x,y] = vpasolve(f3,f2);
x=double(x);
y=double(y);
OcP2=[0,L1,0];
OcpP2=[-x,L1-y,0];
fangxiang=cross(OcP2,OcpP2);

cosr=dot( OcP2,OcpP2)/(norm( OcP2)*norm(OcpP2));
if(fangxiang(3)>=0)
    r=acos(cosr);
else
    r=-acos(cosr);
end
% syms a;
% OcOcp = [x;y];%所有都在Oc系下
% OOc = [sqrt(3)*L0/3+sqrt(3)*L1/6;-L1/2];
% OOcz = [L1/2;sqrt(3)*L0/3+sqrt(3)*L1/6];
% OOcp = OOc+OcOcp;
% OOc=OOc./sqrt(sum(OOc.^2));
% OOcz=OOcz./sqrt(sum(OOcz.^2));
% OOcp=double(OOcp./sqrt(sum(OOcp.^2)));
% f1 = cos(a)*OOc(1)+sin(a)*OOcz(1)-OOcp(1);
% f2 = cos(a)*OOc(2)+sin(a)*OOcz(2)-OOcp(2);
% r1 =double(solve(f1,"real",true));
% r2=double(solve(f2,"real",true));
% %寻找公共解
% rp1=[r1;r2];
% Rp1=sort(rp1);
% dr=diff(Rp1);
% dr=abs(dr);
% [~,ii]=sort(dr);
% r=(Rp1(ii(1))+Rp1(ii(1)+1))/2;
% %最后这个为啥就无解呢啊啊啊啊啊啊啊啊
%把x和y分别带进去求解算的解完全不一样啊
% %%图像
% aa=0:0.01:8;
% yy1=cos(aa).*OOc(1)+sin(aa).*OOcz(1)-OOcp(1);
% yy2=cos(aa).*OOc(2)+sin(aa).*OOcz(2)-OOcp(2);
% plot(aa,yy1)
% hold on
% plot(aa,yy2)
% plot(aa,zeros(length(aa),1))
% hold off
end