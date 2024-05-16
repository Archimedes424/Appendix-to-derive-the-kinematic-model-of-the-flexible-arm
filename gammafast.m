function [r,x,y] = gammafast(h1,h2,h3,L13f,L23f,L0,L1)
%P1P3=(x13,y13)
p13mo=L13f;
p23mo=L23f;
l0=L0;
x13=(2*h1^2*l0^2 - 2*h1*h2*l0^2 - 2*h1*h3*l0^2 + 2*h2*h3*l0^2 + 3*l0^4 + sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*sqrt(l0^4*(4*h1^2 - 4*h1*h2 - 4*h1*h3 + 4*h2^2 - 4*h2*h3 + 4*h3^2 + 9*l0^2)))*sqrt(3)/(4*(h1^2 - 2*h1*h2 + h2^2 + 3*l0^2)*l0);
y13=(2*sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*h1^2 - 2*sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*h1*h2 - 2*sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*h1*h3 + 2*sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*h2*h3 + 3*sqrt(4*h1^2 - 8*h1*h2 + 4*h2^2 + 9*l0^2)*l0^2 - 3*sqrt(l0^4*(4*h1^2 - 4*h1*h2 - 4*h1*h3 + 4*h2^2 - 4*h2*h3 + 4*h3^2 + 9*l0^2)))/(4*(h1^2 - 2*h1*h2 + h2^2 + 3*l0^2));
% x13=max(x1301,x1302);
% y13=min(y1301,y1302);

%//在O3系中用120度关系求解Oc'！r0好像没啥用啊

[x,y] = niu([0,0],L0,L1,x13,y13);
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