function T=exp_ksi(ksi)
theta=sqrt(ksi(4)^2+ksi(5)^2+ksi(6)^2);
rou=[ksi(1);
     ksi(2);
     ksi(3)];
if theta==0
    a=[ksi(4);
       ksi(5);
       ksi(6)];
   t=rou;
else
    a=[double(ksi(4)/theta);
       double(ksi(5)/theta);
       double(ksi(6)/theta)];
    rou=[ksi(1);
         ksi(2);
         ksi(3)];
    J=double(sin(theta))/theta*eye(3,3)+double(1-sin(theta))/theta*(a*a')+double(1-cos(theta))/theta*a';
    t=J*rou;
end
R=v2R(theta,a);
T=[R,t;
    zeros(1,3),1];
end
