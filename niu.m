function [x,y]=niu(chu,L0,L1,x13,y13)
x2=chu(1);
y2=chu(2);
l1=L1;
l0=L0;
f1=1;
f2=1;
P1P3 = [x13,y13];
P1P2 = [sqrt(3)*L0/2,L0/2+L1]; % p.s.全用的横着的！
dx=0.00001;
dy=0.00001;
while(f1*f1+f2*f2>1e-20)
    x1=x2;
    y1=y2;
    OcpP1 = [-sqrt(3)*L0/2-x1,-L0/2-y1];
    f1=norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
    f2=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);
    OcpP1 = [-sqrt(3)*L0/2-x1-dx,-L0/2-y1];
    f1dx=norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
    f2dx=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);
    OcpP1 = [-sqrt(3)*L0/2-x1,-L0/2-y1-dy];
    f1dy=norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
    f2dy=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);
    %     f1=x1*x1-y1-4;
    %     f2=x1-y1*y1-5;
    ykb=[(f1dx-f1)/dx,(f1dy-f1)/dy;(f2dx-f2)/dx,(f2dy-f2)/dy];
    %ykb=[2*x1,-1;1,-2*y1];
    xxj=[x1;y1]-1*ykb\[f1;f2];
    x2z=xxj(1);
    y2z=xxj(2);
    lam=1;
    x2=lam*x2z+(1-lam)*x1;
    y2=lam*y2z+(1-lam)*y1;
    OcpP1 = [-sqrt(3)*L0/2-x2,-L0/2-y2];
    f1x=norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
    f2x=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);
%下山循环，需要时启动
%     while(max(abs([f1x,f2x])-abs([f1,f2]))>0.01)
%         lam=lam/2;
%         lam
%         x2=lam*x2z+(1-lam)*x1;
%         y2=lam*y2z+(1-lam)*y1;
%         OcpP1 = [-sqrt(3)*L0/2-x2,-L0/2-y2];
%         f1x=norm(OcpP1+P1P3)*norm(OcpP1+P1P2)/2+dot(OcpP1+P1P3,OcpP1++P1P2);
%         f2x=norm(OcpP1)*norm(OcpP1+P1P2)/2+dot(OcpP1,OcpP1++P1P2);
%         %         f1x=x2*x2-y2-4;
%         %         f2x=x2-y2*y2-5;
%     end
end
x=x2;
y=y2;
end