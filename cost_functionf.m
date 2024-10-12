function [e,g]=cost_functionf(ksi,A,B,C,U)
e=0.5*ksi'*A*ksi+B*ksi+(ksi)'*U*ksi+0.5*C;
if nargout>1
   g=0.5*(A+A'+4*U)*ksi+B';
end
