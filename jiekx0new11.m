function kx = jiekx0new11(h1,h2,h3,h1t,h2t,h3t,dh1,dh2,dh3,dt,q)
h1d=q(1);
h2d=q(2);
h3d=q(3);
ha0=pqqd11(h1,h2,h3,h1d,h2d,h3d,dh1,dh2,dh3,dt,q);
hah1=pqqd11(h1+dh1,h2,h3,h1d,h2d,h3d,dh1,dh2,dh3,dt,q);
hah2=pqqd11(h1,h2+dh2,h3,h1d,h2d,h3d,dh1,dh2,dh3,dt,q);
hah3=pqqd11(h1,h2,h3+dh3,h1d,h2d,h3d,dh1,dh2,dh3,dt,q);
c1=(hah1-ha0)/dh1;
c2=(hah2-ha0)/dh2;
c3=(hah3-ha0)/dh3;
kx=[c1,c2,c3,zeros(6,6)]*q;
end