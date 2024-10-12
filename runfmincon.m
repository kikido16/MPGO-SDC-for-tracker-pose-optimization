function [ksi_result,fval,history,searchdir] = runfmincon(ksi,A,B,C,U,d1,d2,d3,D1,D2,D3,T2_noise)
 
% Set up shared variables with outfun
history.ksi = [];
history.fval = [];
searchdir = [];
t2=[-100;0;0];
angle2=[0;0;0];
R2=rpy2r(angle2(1),angle2(2),angle2(3),'degree');
trans_errors=[];
rot_errors=[];
% Call optimization

options = optimoptions('fmincon','ConstraintTolerance',0.05,'FunctionTolerance',1e-3,...
        'MaxFunctionEvaluations',5000,'MaxIterations',1000,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(ksi,lambda)hessianfcn(ksi,lambda,A,D1,D2,D3,U),...
         'OutputFcn',@outfun);
[ksi_result,fval,exitflag,output]=fmincon(@(ksi)objfun(ksi,A,B,C,U),ksi,[],[],[],[],[],[],@(ksi)confun(ksi,d1,d2,d3,D1,D2,D3),options);
 
 function stop = outfun(ksi,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.ksi = [history.ksi; ksi];
           it_times=length(history.fval);
               ksi_itr=history.ksi(18*(it_times-1)+6:18*(it_times-1)+11);
               T2_itr=SE3.exp(ksi_itr)*T2_noise;
               t2_itr=T2_itr.t
               R2_itr=T2_itr.SO3.R;
               trans_error=norm(t2-t2_itr);
               rot_error=acos(double((trace(R2*R2_itr')-1))/2);
               trans_errors(it_times)=trans_error;
               rot_errors(it_times)=rot_error;
         % Concatenate current search direction with 
         % searchdir.
%            searchdir = [searchdir;... 
%                         optimValues.searchdirection'];
           plot(trans_error,'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'.
%            text(ksi(1)+.15,ksi(2),... 
%                 num2str(optimValues.iteration));
           title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
 end
 
function [e,g]=objfun(ksi,A,B,C,U)
e=0.5*ksi'*A*ksi+B*ksi+(ksi)'*U*ksi+0.5*C;
if nargout>1
   g=0.5*(A+A'+4*U)*ksi+B';
end
end
 
function [res,ceq,gc,geq]=confun(ksi,d1,d2,d3,D1,D2,D3)
res1=ksi'*(D1)'*D1*ksi+2*d1'*D1*ksi+d1'*d1;
res2=ksi'*(D2)'*D2*ksi+2*d2'*D2*ksi+d2'*d2;
res3=ksi'*(D3)'*D3*ksi+2*d3'*D3*ksi+d3'*d3;
res=[res1-10002;
    -res1+9998;
    res2-20003;
    -res2+19997;
    res3-10002;
    -res3+9998]; 
ceq=[];
if nargout>2
    r1g=((D1)'*D1+((D1)'*D1)')*ksi+2*(d1'*D1)';
    r2g=((D2)'*D2+((D2)'*D2)')*ksi+2*(d2'*D2)';
    r3g=((D3)'*D3+((D3)'*D3)')*ksi+2*(d3'*D3)';
    gc=[r1g,-r1g,r2g,-r2g,r3g,-r3g];
    geq=[];
end
end
end