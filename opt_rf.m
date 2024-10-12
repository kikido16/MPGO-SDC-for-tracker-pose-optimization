function T=opt_rf(Ti_noise,Tij_noise,sigma)
%init values
T1_noise=Ti_noise(1);
T2_noise=Ti_noise(2);
T3_noise=Ti_noise(3);
T4_noise=Ti_noise(4);
T12_noise=Tij_noise(1);
T23_noise=Tij_noise(2);
T24_noise=Tij_noise(3);
T34_noise=Tij_noise(4);
T13_noise=Tij_noise(5);
T14_noise=Tij_noise(6);
delta_ksi=1;
it_times=0;

u=1;
u_trans=1;
u_rot=1;
u_ksi=[u_trans,u_trans,u_trans,u_rot,u_rot,u_rot];

u1=20;
v=2;
fval_current=0;
while (delta_ksi>5e-3)
    
    if (it_times>40)
        break;
    end
    uu=u*[u_ksi,u_ksi,u_ksi];
    U=diag(uu);
    it_times=it_times+1;
    e12=T12_noise/(T2_noise)*T1_noise;
    e23=T23_noise/(T3_noise)*T2_noise;
    e24=T24_noise/(T4_noise)*T2_noise;
    e34=T34_noise/(T4_noise)*T3_noise;
    e13=T13_noise/(T3_noise)*T1_noise;
    e14=T14_noise/(T4_noise)*T1_noise;
    I=eye(6);
    eta12=e12.logs;
    eta23=e23.logs;
    eta24=e24.logs;
    eta34=e34.logs;
    eta13=e13.logs;
    eta14=e14.logs;
    e=[norm(eta12),norm(eta23),norm(eta24),norm(eta34),norm(eta13),norm(eta14)];
    
    %calculate robust weights
    w=welsch_weight(e,u1);
    
    % cost function linearization
    rou12=[eta12(1),eta12(2),eta12(3)];
    phi12=[eta12(4),eta12(5),eta12(6)];
    rou23=[eta23(1),eta23(2),eta23(3)];
    phi23=[eta23(4),eta23(5),eta23(6)];
    rou24=[eta24(1),eta24(2),eta24(3)];
    phi24=[eta24(4),eta24(5),eta24(6)];
    
    rou34=[eta34(1),eta34(2),eta34(3)];
    phi34=[eta34(4),eta34(5),eta34(6)];
    rou13=[eta13(1),eta13(2),eta13(3)];
    phi13=[eta13(4),eta13(5),eta13(6)];
    rou14=[eta14(1),eta14(2),eta14(3)];
    phi14=[eta14(4),eta14(5),eta14(6)];

    Jri12=I+0.5*[skew(phi12),skew(rou12);zeros(3,3),skew(phi12)];
    Jri23=I+0.5*[skew(phi23),skew(rou23);zeros(3,3),skew(phi23)];
    Jri24=I+0.5*[skew(phi24),skew(rou24);zeros(3,3),skew(phi24)];
    Jri34=I+0.5*[skew(phi34),skew(rou34);zeros(3,3),skew(phi34)];
    Jri13=I+0.5*[skew(phi13),skew(rou13);zeros(3,3),skew(phi13)];
    Jri14=I+0.5*[skew(phi14),skew(rou14);zeros(3,3),skew(phi14)];
    iT1_noise=inv(T1_noise);
    iT2_noise=inv(T2_noise);
    iT3_noise=inv(T3_noise);
    iT4_noise=inv(T4_noise);

    J12=[-Jri12*iT1_noise.Ad(),zeros(6,6),zeros(6,6)];
    J23=[Jri23*iT2_noise.Ad(),-Jri23*iT2_noise.Ad(),zeros(6,6)];
    J24=[Jri24*iT2_noise.Ad(),zeros(6,6),-Jri24*iT2_noise.Ad()];
    
    J34=[zeros(6,6),Jri34*iT3_noise.Ad(),-Jri34*iT3_noise.Ad()];
    J13=[zeros(6,6),-Jri13*iT1_noise.Ad(),zeros(6,6)];
    J14=[zeros(6,6),zeros(6,6),-Jri14*iT1_noise.Ad()];
    
    J=[Jri12*iT1_noise.ad(),-Jri12*iT1_noise.ad(),Jri23*iT2_noise.ad(),-Jri23*iT2_noise.ad(),Jri24*iT2_noise.ad(),-Jri24*iT2_noise.ad()];
    A12=J12'*J12;
    A23=J23'*J23;
    A24=J24'*J24;
    A34=J34'*J34;
    A13=J13'*J13;
    A14=J14'*J14;
    
    B12=eta12*J12;
    B23=eta23*J23;
    B24=eta24*J24;
    B34=eta34*J34;
    B13=eta13*J13;
    B14=eta14*J14;
    
    C12=eta12*eta12';
    C23=eta23*eta23';
    C24=eta24*eta24';
    C34=eta34*eta34';
    C13=eta13*eta13';
    C14=eta14*eta14';
    
    A=w(1)*A12+w(2)*A23+w(3)*A24+w(4)*A34+w(5)*A13+w(6)*A14;
    B=w(1)*B12+w(2)*B23+w(3)*B24+w(4)*B34+w(5)*B13+w(6)*B14;
    C=w(1)*C12+w(2)*C23+w(3)*C24+w(4)*C34+w(5)*C13+w(6)*C14;
    if (it_times==1)
        u=1e-4*max(diag(A));
    end
    % restrciton linearization
    B1=[0;0;100];
    B2=[0;0;100];
    B3=[0;0;100];
    P1=iT1_noise*B1;
    P2=iT3_noise*B2;
    P3=iT4_noise*B3;
    r1=norm(P1-P2);
    r2=norm(P1-P3);
    r3=norm(P2-P3);
    dres1=abs(1e4-(r1)^2+2e4-(r2)^2+1e4-(r3)^2);
    DO1=iT1_noise.T()*[eye(3),-skew(B1);zeros(1,6)];
    DO2=iT3_noise.T()*[eye(3),-skew(B2);zeros(1,6)];
    DO3=iT4_noise.T()*[eye(3),-skew(B3);zeros(1,6)];
    d1=[P1-P2;0];
    d2=[P1-P3;0];
    d3=[P2-P3;0];
    D1=[zeros(4,6),DO2,zeros(4,6)];
    D2=[zeros(4,6),zeros(4,6),DO3];
    D3=[zeros(4,6),-DO2,DO3];
    ksi=zeros(18,1);
 
    %optimization settings
    options = optimoptions('fmincon','ConstraintTolerance',0.05,'FunctionTolerance',1e-3,...
        'MaxFunctionEvaluations',5000,'MaxIterations',1000,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(ksi,lambda)hessianfcn(ksi,lambda,A,D1,D2,D3,U));
    options.CheckGradients=false;
    [ksi_result,fval,exitflag,output]=fmincon(@(ksi)cost_functionf(ksi,A,B,C,U),ksi,[],[],[],[],[],[],@(ksi)ball_resf(ksi,d1,d2,d3,D1,D2,D3),options);
    error_plot(it_times)=fval;
    restrictions_result=ball_resf(ksi_result,d1,d2,d3,D1,D2,D3);    
    delta_ksi=norm(ksi_result);
    
    %calculate the possible optimized pose
    ksi2=ksi_result(1:6);
    ksi3=ksi_result(7:12);
    ksi4=ksi_result(13:18);
    T2_new=SE3.exp(ksi2)*T2_noise;
    T3_new=SE3.exp(ksi3)*T3_noise;
    T4_new=SE3.exp(ksi4)*T4_noise;
    P1_new=iT1_noise*B1;
    P2_new=inv(T3_new)*B2;
    P3_new=inv(T4_new)*B3;
    r1_new=norm(P1_new-P2_new);
    r2_new=norm(P1_new-P3_new);
    r3_new=norm(P2_new-P3_new);

    %paras update
    delta_ek=fval-fval_current;
    fval_current=fval;
    grad_delta_ek=(0.5*(A+A'+4*U)*ksi_result+B')'*ksi_result;
    if(it_times==1)
        rou=abs(delta_ek/grad_delta_ek);
    else
        rou=delta_ek/grad_delta_ek;
    end
    if (rou>0&&u>0.1)
        u=u*max(1/3,1-(2*rou-1)^3);
        v=2;
        T2_noise=SE3.exp(ksi2)*T2_noise;
        T3_noise=SE3.exp(ksi3)*T3_noise;
        T4_noise=SE3.exp(ksi4)*T4_noise;
    else
        u=u*v;
        v=2*v;
    end

    
    if(mod(it_times,4)==0)
        u1 = u1/2;
    end
end
T=[T1_noise,T2_noise,T3_noise,T4_noise];
end