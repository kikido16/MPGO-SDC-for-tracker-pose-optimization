function Hout=hessianfcn(ksi,lambda,A,D1,D2,D3,U)
H=0.5*(A+A'+4*U);
Hg1=(D1)'*D1+((D1)'*D1)';
Hg2=(D2)'*D2+((D2)'*D2)';
Hg3=(D3)'*D3+((D3)'*D3)';

Hout=H+lambda.ineqnonlin(1)*Hg1+lambda.ineqnonlin(2)*(-Hg1)+...
    lambda.ineqnonlin(3)*Hg2+lambda.ineqnonlin(4)*(-Hg2)+...
    lambda.ineqnonlin(5)*Hg3+lambda.ineqnonlin(6)*(-Hg3);
end
