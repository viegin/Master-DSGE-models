// new keynesian model with inertial and taylor rule 
// Master Macroeconomic 2020

//-------------------------------------------------
// Variables
//-------------------------------------------------

var p, y, i, r, g, u, v;
varexo eg eu ev;

//-------------------------------------------------
// Parameters
//-------------------------------------------------

parameters beta,  sigma, omega, alpha, kappa, h ,gamma, d, rhog, rhou, rhov;

//-------------------------------------------------
// Parameter Values
//-------------------------------------------------

beta = 0.99;
sigma = 1;
omega = 0.5;
alpha = 3;
kappa = (1-omega)*(1-beta*omega)/(alpha*omega);
h = 0.5;
gamma = 0.5; 
d =1.5;
rhog = 0.5; 
rhou = 0.5; 
rhov = 0.5;

//-------------------------------------------------
// Model Equations
//-------------------------------------------------

model;
i = d*p + v;
y = h*y(-1)+(1-h)*y(+1) - (1/sigma)*(i - p(+1))+ g;
p = gamma*p(-1)+(1-gamma)*p(+1) + kappa*y + u; 
r = i - p(+1);

v = rhov*v(-1)+ ev;
u = rhou*u(-1)+ eu;
g = rhog*g(-1)+ eg;

end;

//-------------------------------------------------
// Initial Values and Steady State
//-------------------------------------------------

initval;
p = 0;
y = 0;
u = 0;
end;

//-------------------------------------------------
// Shocks - stadard errors
//-------------------------------------------------

shocks;
var ev; stderr 1;
var eu; stderr 1;
var eg; stderr 1;
end;

//-------------------------------------------------
//  Start stochastic simulation
//-------------------------------------------------


stoch_simul(linear, irf=20)i y p r;

