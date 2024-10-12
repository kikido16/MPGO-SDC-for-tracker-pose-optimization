function [res,ceq,gc,geq]=ball_resf(ksi,d1,d2,d3,D1,D2,D3)
% ksi1=ksi(1:6);
% ksi1=zeros(6,1);
% ksi3=ksi(7:12);
% ksi4=ksi(13:18);
% C1=[P1-P2;0];
% C2=[P1-P3;0];
% D1=norm(C1+DO2*ksi3);
% D2=norm(C2+DO3*ksi4);
% D1=norm((C1-(DO1*ksi1-DO2*ksi3)));
% D2=norm((C2-(DO1*ksi1-DO3*ksi4)));
res1=ksi'*(D1)'*D1*ksi+2*d1'*D1*ksi+d1'*d1;
res2=ksi'*(D2)'*D2*ksi+2*d2'*D2*ksi+d2'*d2;
res3=ksi'*(D3)'*D3*ksi+2*d3'*D3*ksi+d3'*d3;

res=[res1-10002;
    -res1+9998;
    res2-20003;
    -res2+19997;
    res3-10002;
    -res3+9998]; 
% res=[(D1)-99.99;
%     -(D1)+99.99;
%     (D2)-141.41;
%     -(D2)+141.41];
% res=[ksi1'*(DO1)'*DO1*ksi1+ksi3'*(DO2)'*DO2*ksi3-2*C1'*DO1*ksi1+2*C1'*DO1*ksi3-2*(DO1*ksi1)'*DO2*ksi3+C1'*C1-10002;
%     10002-(ksi1'*(DO1)'*DO1*ksi1+ksi3'*(DO2)'*DO2*ksi3-2*C1'*DO1*ksi1+2*C1'*DO1*ksi3-2*(DO1*ksi1)'*DO2*ksi3+C1'*C1);
%     ksi1'*(DO1)'*DO1*ksi1+ksi4'*(DO3)'*DO3*ksi4-2*C2'*DO1*ksi1+2*C2'*DO1*ksi4-2*(DO1*ksi1)'*DO3*ksi4+C2'*C2-20003;
%     20003-(ksi1'*(DO1)'*DO1*ksi1+ksi4'*(DO3)'*DO3*ksi4-2*C2'*DO1*ksi1+2*C2'*DO1*ksi4-2*(DO1*ksi1)'*DO3*ksi4+C2'*C2)];
ceq=[];
if nargout>2
    r1g=((D1)'*D1+((D1)'*D1)')*ksi+2*(d1'*D1)';
    r2g=((D2)'*D2+((D2)'*D2)')*ksi+2*(d2'*D2)';
    r3g=((D3)'*D3+((D3)'*D3)')*ksi+2*(d3'*D3)';
    gc=[r1g,-r1g,r2g,-r2g,r3g,-r3g];
%     gc=[r1g,-r1g,r2g,-r2g];
    geq=[];
end
% if nargout>4
%     H1=(D1)'*D1+((D1)'*D1)';
%     H2=(D2)'*D2+((D2)'*D2)';
%     Hc=[H1,-H1,H2,-H2];
%     Heq=[];
% end
