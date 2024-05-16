function [h1m,h2m,h3m]=x2hm(xm,ym,zm)
hm0=[zm;zm;zm];
hm1=hm0;
om=[xm;ym;zm];
h0 = 3;
L0 = 2;
%设置关键参数
ex=100;
ej=100;
es=[];
r=0.0001;
while(ex>r)
    lam=1;
    ykb=suanykb(hm1(1),hm1(2),hm1(3),h0,L0);
    yhm1=jieA(hm1(1),hm1(2),hm1(3),h0,L0)*[0;0;0;1];
    yhm1=yhm1(1:3)-om;
    hm2=hm1-1*ykb\yhm1;
    ygu=jieA(hm2(1),hm2(2),hm2(3),h0,L0)*[0;0;0;1];
    ygu=ygu(1:3)-om;
    ex=norm(ygu);
    while(ex>ej)
        lam=lam/2;
        hm2s=hm2*lam+hm1*(1-lam);
        ygu=jieA(hm2s(1),hm2s(2),hm2s(3),h0,L0)*[0;0;0;1];
        ygu=ygu(1:3)-om;
        ex=norm(ygu);
    end
    hm2=hm2s;
    es=[es;ex];
    hm1=hm2;
    ej=ex;
end
h1m=hm2(1);
h2m=hm2(2);
h3m=hm2(3);
end