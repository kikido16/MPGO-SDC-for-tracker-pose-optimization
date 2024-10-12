clear all;
% initial values
init_trans_error=[];
init_rot_error=[];
trans_noise_level=0.1;
rot_noise_level=0.005;

% different noise levels
for r_noise=1:5
    for t_noise=1:5
        for i=1:20
            angle1=[0;0;0];
            angle2=[0;0;0];
            angle3=[0;0;0];
            angle4=[0;0;0];
            r12=[0;0;0];
            r23=[0;0;0];
            r24=[0;0;0];
            R1=rpy2r(angle1(1),angle1(2),angle1(3),'degree');
            R2=rpy2r(angle2(1),angle2(2),angle2(3),'degree');
            R3=rpy2r(angle3(1),angle3(2),angle3(3),'degree');
            R4=rpy2r(angle4(1),angle4(2),angle4(3),'degree');
            R12=rpy2r(r12(1),r12(2),r12(3),'degree');
            R23=rpy2r(r23(1),r23(2),r23(3),'degree');
            R24=rpy2r(r24(1),r24(2),r24(3),'degree');
            
            %set initial pose values
            t1=[0;0;-100];
            t2=[-100;0;0];
            t3=[-100;0;-100];
            t4=[-100;0;-200];
            t12=[-100;0;100];
            t23=[0;0;-100];
            t24=[0;0;-200];
            T1=SE3(R1,t1);
            T2=SE3(R2,t2);
            T3=SE3(R3,t3);
            T4=SE3(R4,t4);
            T12=SE3(R12,t12);
            T23=SE3(R23,t23);
            T24=SE3(R24,t24);
            
            %set random noise values under the noise level
            sigma=[trans_noise_level,trans_noise_level,trans_noise_level,...
                rot_noise_level,rot_noise_level,rot_noise_level];          
            rng(i);
            xi1=sigma.*randn(1,6);
            xi2=sigma.*randn(1,6);
            xi3=sigma.*randn(1,6);
            xi4=sigma.*randn(1,6);
            xi5=sigma.*randn(1,6);
            xi6=sigma.*randn(1,6);
            
            % calculate perturbed relative poses
            T1_noise=T1;
            T12_noise=SE3.exp(xi1)*T12;
            T23_noise=SE3.exp(xi2)*T23;
            T24_noise=SE3.exp(xi3)*T24;
            T34_noise=inv(T23_noise)*T24_noise;
            T13_noise=T12_noise*T23_noise;
            T14_noise=T12_noise*T24_noise;
            
            % calculate perturbed absolute poses
            T2_noise=T1_noise*T12_noise;
            T3_noise=T2_noise*T23_noise;
            T4_noise=T2_noise*T24_noise;
            R2_noise=T2_noise.SO3.R;
            R3_noise=T3_noise.SO3.R;
            R4_noise=T4_noise.SO3.R;
            t2_noise=T2_noise.t;
            t3_noise=T3_noise.t;
            t4_noise=T4_noise.t;
            
            %calculate init tracker pose errors
            init_rot_error2=acos(double((trace(R2*R2_noise')-1))/2);
            init_trans_error2=norm(t2-t2_noise);
            init_trans_error(i)=init_trans_error2;
            init_rot_error(i)=init_rot_error2;

            Ti_noise=[T1_noise,T2_noise,T3_noise,T4_noise];
            Tij_noise=[T12_noise,T23_noise,T24_noise,T34_noise,T13_noise,T14_noise];
           
            %optimization with robust function
            opt_rf_resultT=opt_rf(Ti_noise,Tij_noise,sigma);
            [opt_rf_trans_error,opt_rf_rot_error]=err_cal(opt_rf_resultT);
            opt_rf_trans_errors(i)=opt_rf_trans_error(2);
            opt_rf_rot_errors(i)=opt_rf_rot_error(2);
                                            
            % %optimization with covariance
            opt_cov_resultT=opt_cov(Ti_noise,Tij_noise,sigma);
            [opt_cov_trans_error,opt_cov_rot_error]=err_cal(opt_cov_resultT);
            opt_cov_trans_errors(i)=opt_cov_trans_error(2);
            opt_cov_rot_errors(i)=opt_cov_rot_error(2);
            
        end
        %trans and rot error calculation
        avg_init_trans_error=mean(init_trans_error);
        avg_opt_cov_trans_error=mean(opt_cov_trans_errors);
        avg_opt_rf_trans_error=mean(opt_rf_trans_errors);
        
        avg_init_rot_error=mean(init_rot_error);
        avg_opt_cov_rot_error=mean(opt_cov_rot_errors);
        avg_opt_rf_rot_error=mean(opt_rf_rot_errors);
                
        avg_init_trans_errors(r_noise,t_noise)=avg_init_trans_error;
        avg_init_rot_errors(r_noise,t_noise)=avg_init_rot_error;
        cov_trans_errors(r_noise,t_noise)=avg_opt_cov_trans_error;
        rf_trans_errors(r_noise,t_noise)=avg_opt_rf_trans_error;
        cov_rot_errors(r_noise,t_noise)=avg_opt_cov_rot_error;
        rf_rot_errors(r_noise,t_noise)=avg_opt_rf_rot_error;

        
        trans_noise_level=trans_noise_level+0.05;
        rot_noise_level=rot_noise_level+0.005;
    end
end

