function b = cosbeta1hanshu(h1,h2,h3,L0)
p1p2=[sqrt(3)*L0/2,3*L0/2,h2-h1];
 xita1=[0,sqrt(3)*(h2-h1)*L0,-3*sqrt(3)*L0*L0/2];

 xita2=[(3*(h3-h1)*L0)/2,sqrt(3)*((h2-h1)*L0-((h3-h1)*L0)/2),-3*sqrt(3)*L0*L0/2];
 aa=cross(xita1,xita2);
 cb=dot(xita1,xita2)/(norm(xita1)*norm(xita2));
 if(dot(p1p2,aa)>=0)
     b=acos(cb);
 else
     b=-acos(cb);
 end
 b=real(b);
% n = 3*(h2-h1)*(h2-h1)-3*(h3-h1)*(h2-h1)/2+27*L0*L0/4;
%  m = sqrt(3*(h2-h1)*(h2-h1)+27*L0*L0/4)*sqrt(9*(h2-h1)*(h2-h1)/4+27*L0*L0/4+1.732*(h2-h3/2-h1/2));
%  b = acos(n/m)-pi/2;

end
