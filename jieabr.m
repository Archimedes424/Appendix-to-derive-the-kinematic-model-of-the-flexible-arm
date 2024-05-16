function [alpha,b,r] = jieabr(h1,h2,h3,L0)
alpha = atan((h2-h1)/(3*L0/2));
b = cosbeta1hanshu(h1,h2,h3,L0);
L1 = sqrt(9*L0*L0/4+(h2-h1)*(h2-h1))-L0/2;
L13f = 3*L0*L0+(h3-h1)*(h3-h1);
L23f = 3*L0*L0+(h3-h2)*(h3-h2);  %L13f即L13的平方
[r,~,~] = gammafast(h1,h2,h3,L13f,L23f,L0,L1);
end