// new keynesian model with taylor rule

var i, x, p;
varexo v g u;

parameters psi, beta, lambda, alpha, a, b, rhov, rhog, rhou;

psi = 0.1;
beta = 0.99;
lambda = 0.2;
alpha = 0.85;
a = 0.5 ;
b = 1.5;
rhov = 0.01;
rhog = 0.02;
rhou = 0.02;

model(linear);
i = a*x + b*p + v;
x = alpha*x(-1)+(1-alpha)*x(+1) - psi*(i - p(+1))+ g;
p = beta*(alpha*p(-1)+(1-alpha)*p(+1)) + lambda*x + u;

end;

initval;
i = 0;
x = 0;
p = 0;
end;

steady;

shocks;
var v; stderr rhov;
var g; stderr rhog;
var u; stderr rhou;
end;

stoch_simul(order=1,drop=0,periods=100,irf=10,nomoments) x p;

save data1 x p;

estimated_params;
  a,normal_pdf,8,0.5;
  b,normal_pdf,1.5,0.5;
  alpha, gamma_pdf, 0.8, 0.3;
  stderr v,gamma_pdf,0.02,0.015;
  stderr g,gamma_pdf,0.02,0.015;
  stderr u,gamma_pdf,0.02,0.015;
 end;


varobs x p;

estimation(datafile=data1 ,mode_compute=4,mh_nblocks=1,mode_check,mh_jscale=1,mh_replic=3000);
