function w = welsch_weight(r,u)

w = exp(-r.^2./(2*u^2));

end