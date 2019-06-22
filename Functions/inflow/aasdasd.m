
rmean = 0.5;
phi1 = sqrt(3); phi2 = sqrt(7); phi3 = sqrt(11); phi4 = sqrt(30)/2; phi5 = sqrt(5);
phi6 = sqrt(210)/4; phi7 = sqrt(1155/8); phi8 = 3/4*sqrt(35); phi9 = sqrt(3465/128);
a(1) = phi1;
a(2) = phi2*(1-5/2*rmean^2);
a(3) = phi3*(1-7*rmean^2+63/8*rmean^4);
a(4) = phi4*rmean;
a(5) = phi5*(3*rmean-21/4*rmean^3);
a(6) = phi6*rmean^2;
a(7) = phi7*(rmean^2-1.5*rmean^4);
a(8) = phi8*rmean^3;
a(9) = phi9*rmean^4;
a(10) = phi4*rmean;
a(11) = phi5*(3*rmean-21/4*rmean^3);
a(12) = phi6*rmean^2;
a(13) = phi7*(rmean^2-1.5*rmean^4);
a(14) = phi8*rmean^3;
a(15) = phi9*rmean^4;
            
ali(1) = RadialShapeFunc(0,1, rmean);
ali(2) = RadialShapeFunc(0,3, rmean);
ali(3) = RadialShapeFunc(0,5, rmean);

ali(4) = RadialShapeFunc(1,2, rmean);
ali(5) = RadialShapeFunc(1,4, rmean);
ali(6) = RadialShapeFunc(2,3, rmean);
ali(7) = RadialShapeFunc(2,5, rmean);
ali(8) = RadialShapeFunc(3,4, rmean);
ali(9) = RadialShapeFunc(4,5, rmean);

ali(10) = RadialShapeFunc(1,2, rmean);
ali(11) = RadialShapeFunc(1,4, rmean);
ali(12) = RadialShapeFunc(2,3, rmean);
ali(13) = RadialShapeFunc(2,5, rmean);
ali(14) = RadialShapeFunc(3,4, rmean);
ali(15) = RadialShapeFunc(4,5, rmean);


