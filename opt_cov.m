function T=opt_cov(Ti_noise,Tij_noise,sigma)
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

B1=[0;0;100];
B2=[0;0;100];
B3=[0;0;100];
u=5;
u_trans=1;
u_rot=1;
u_ksi=[u_trans,u_trans,u_trans,u_rot,u_rot,u_rot];

v=2;
fval_current=0;
    while (delta_ksi>5e-3)
        % cost function linearization
        it_times=it_times+1;
        if (it_times>40)
            break;
        end
        uu=u*[u_ksi,u_ksi,u_ksi];
        U=diag(uu);
        cov=(diag(abs(sigma)))^2;
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

        J12=sqrt(inv(cov))*[-Jri12*iT1_noise.Ad(),zeros(6,6),zeros(6,6)];
        J23=sqrt(inv(cov))*[Jri23*iT2_noise.Ad(),-Jri23*iT2_noise.Ad(),zeros(6,6)];
        J24=sqrt(inv(cov))*[Jri24*iT2_noise.Ad(),zeros(6,6),-Jri24*iT2_noise.Ad()];

        J34=inv(cov)*[zeros(6,6),Jri34*iT3_noise.Ad(),-Jri34*iT3_noise.Ad()];
        J13=inv(cov)*[zeros(6,6),-Jri13*iT1_noise.Ad(),zeros(6,6)];
        J14=inv(cov)*[zeros(6,6),zeros(6,6),-Jri14*iT1_noise.Ad()];
        J=[Jri12*iT1_noise.ad(),-Jri12*iT1_noise.ad(),Jri23*iT2_noise.ad(),-Jri23*iT2_noise.ad(),Jri24*iT2_noise.ad(),-Jri24*iT2_noise.ad()];

        A12=J12'*J12;
        A23=J23'*J23;
        A24=J24'*J24;
        A34=J34'*J34;
        A13=J13'*J13;
        A14=J14'*J14;

        B12=eta12*sqrt(inv(cov))*J12;
        B23=eta23*sqrt(inv(cov))*J23;
        B24=eta24*sqrt(inv(cov))*J24;
        B34=eta34*inv(cov)*J34;
        B13=eta13*inv(cov)*J13;
        B14=eta14*inv(cov)*J14;
        
        C12=eta12*eta12';
        C23=eta23*eta23';
        C24=eta24*eta24';
        C34=eta34*eta34';
        C13=eta13*eta13';
        C14=eta14*eta14';
       
        A=A12+A23+A24+A34+A13+A14;
        B=B12+B23+B24+B34+B13+B14;
        C=C12+C23+C24+C34+C13+C14;
        if (it_times==1)
            u=1e-4*max(diag(A));
        end
        
        % restrciton linearization
        P1=iT1_noise*B1;
        P2=iT3_noise*B2;
        P3=iT4_noise*B3;
        DO1=iT1_noise.T()*[eye(3),-skew(B1);zeros(1,6)];
        DO2=iT3_noise.T()*[eye(3),-skew(B2);zeros(1,6)];
        DO3=iT4_noise.T()*[eye(3),-skew(B3);zeros(1,6)];
        d1=[P1-P2;0];
        d2=[P1-P3;0];
        d3=[P2-P3;0];
        D1=[zeros(4,6),DO2,zeros(4,6)];
        D2=[zeros(4,6),zeros(4,6),DO3];
        D3=[zeros(4,6),-DO2,DO3];
        
        %optimization settings
        ksi=zeros(18,1);
        options = optimoptions('fmincon','ConstraintTolerance',0.05,'FunctionTolerance',1e-3,...
        'MaxFunctionEvaluations',5000,'MaxIterations',1000,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(ksi,lambda)hessianfcn(ksi,lambda,A,D1,D2,D3,U));
        options.CheckGradients=false;
        [ksi_result,fval,exitflag,output]=fmincon(@(ksi)cost_functionf(ksi,A,B,C,U),ksi,[],[],[],[],[],[],@(ksi)ball_resf(ksi,d1,d2,d3,D1,D2,D3),options);
        delta_ksi=norm(ksi_result);
        ksi2=ksi_result(1:6);
        ksi3=ksi_result(7:12);
        ksi4=ksi_result(13:18);
        
        %calculate the possible optimized pose
        T2_new=SE3.exp(ksi2)*T2_noise;
        T3_new=SE3.exp(ksi3)*T3_noise;
        T4_new=SE3.exp(ksi4)*T4_noise;
        e12_new=T12_noise/(T2_new)*T1_noise;
        e23_new=T23_noise/(T3_new)*T2_new;
        e24_new=T24_noise/(T4_new)*T2_new;
        e34_new=T34_noise/(T4_new)*T3_new;
        e13_new=T13_noise/(T3_new)*T1_noise;
        e14_new=T14_noise/(T4_new)*T1_noise;
        eta12_new=e12_new.logs;
        eta23_new=e23_new.logs;
        eta24_new=e24_new.logs;
        eta34_new=e34_new.logs;
        eta13_new=e13_new.logs;
        eta14_new=e14_new.logs;
        
    %paras update
    delta_ek=fval-fval_current;
    fval_current=fval;
    grad_delta_ek=(0.5*(A+A'+4*U)*ksi_result+B')'*ksi_result;
        if(it_times==1)
            rou=abs(delta_ek/grad_delta_ek);
        else
            rou=delta_ek/grad_delta_ek;
        end
        if (rou>0)
            u=u*max(1/3,1-(2*rou-1)^3);
            v=2;
            T2_noise=SE3.exp(ksi2)*T2_noise;
            T3_noise=SE3.exp(ksi3)*T3_noise;
            T4_noise=SE3.exp(ksi4)*T4_noise;
        else
            u=u*v;
            v=2*v;
        end
    end
T=[T1_noise,T2_noise,T3_noise,T4_noise];
end