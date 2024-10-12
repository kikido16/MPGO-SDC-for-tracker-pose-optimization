function J = computeJocabianMatrix(Ksai,m)

    if m==-1
        Ksai=-Ksai;
    end

    w_Ksai=Ksai(4:6);
    v_Ksai=Ksai(1:3);
    theta=norm(w_Ksai,2);
%     adjoint_ou=Ksai.ad();
    w_Ksai_h=slove_antisymmetricmatrix(w_Ksai);
    v_Ksai_h=slove_antisymmetricmatrix(v_Ksai);

    adjoint_ou=[w_Ksai_h v_Ksai_h;zeros(3,3) w_Ksai_h];
    

    coeff_one=(4-theta*sin(theta)-4*cos(theta))/2/theta^2;
    coeff_two=(4*theta-5*sin(theta)+theta*cos(theta))/2/theta^3;
    coeff_three=(2-theta*sin(theta)-2*cos(theta))/2/theta^4;
    coeff_four=(2*theta-3*sin(theta)+theta*cos(theta))/2/theta^5;

    J=eye(6,6)+coeff_one*adjoint_ou+coeff_two*adjoint_ou^2+coeff_three*adjoint_ou^3+coeff_four*adjoint_ou^4;
    
    if theta==0
    
        J=eye(6,6);

    end


end
