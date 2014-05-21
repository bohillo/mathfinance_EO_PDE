1;
function [x,V1,V2,F,S0,dx,dt,DF_d,DF_f] = Parisian_out(F_bid, F_ask, barrier,day_hat, strike,
			   issue_date,expire_date,PPO,OSO,price_type, \
                           barrier_type, payoff_type,isAsian,indsig)
# set parameter
global Mx; 
global Mt;
global FC_DOM;
global DCC;
global dsigma;

if (price_type == "bid")
    F = F_bid;
else
    F = F_ask;
endif

[T, T_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, price_type);

  t = linspace(0,T,Mt+1);
  Tau = year_frac(issue_date, day_shift2(issue_date,FC_DOM,day_hat), DCC);
  r0 = -log(DF_d)/T;
  r1 = -log(DF_f)/T;
  S0 = F * DF_d / DF_f;
  xmax = log(set_xmax(S0, T, vol(0, S0, strike, T)));
  xmin = -xmax;
  x = linspace(xmin, xmax, Mx + 2)';
  dt = T / Mt;
  sigma = vol(t, x, strike, T)+dsigma*indsig;
  bar = log(barrier);
  dx = 2*xmax/(Mx+1);
  NoTau = floor(Tau/dt);


alphap = -sigma(2:Mx,Mt+1).^2/2/dx^2 +(sigma(2:Mx,Mt+1).^2/2-(r0-r1))/2/dx;
alpham = -sigma(2:Mx,Mt+1).^2/2/dx^2 -(sigma(2:Mx,Mt+1).^2/2-(r0-r1))/2/dx;
beta = sigma(2:Mx,Mt+1).^2/dx^2+r0;
A = spdiags([alpham, beta, alphap], -1:1, Mx, Mx);
theta = 0.5;
B = speye(Mx,Mx) + theta*dt*A;
C = speye(Mx,Mx) - (1-theta)*dt*A;

inx = find(x>bar,1,'first');

if strcmp(payoff_type,"call")
  left_cond = zeros(1,length(t));
  u = max(exp(x)-strike,0);
  u = repmat(u,1,NoTau+1);
  u(1,:) = left_cond(1);
  u(:,1) = 0; 
elseif strcmp(payoff_type,"put")
  left_cond = strike*exp(-r0*(T-t));
  u = max(strike-exp(x),0);
  u = repmat(u,1,NoTau+1);	   
  u(1,:) = left_cond(1);
  u(:,1) = 0; 
  
endif

f = zeros(Mx,1);


for m = 1:Mt

  lastu = u;
  j = 2:NoTau+1;
  if strcmp(barrier_type,"up")
	f = C*[repmat(lastu(2:inx-2,NoTau+1),1,NoTau);lastu(inx-1:end-1,j-1)]; 
  elseif strcmp(barrier_type,"down")
	f = C*[lastu(2:inx-2,j-1);repmat(lastu(inx-1:end-1,NoTau+1),1,NoTau)];
  endif
  u(:,j) = zeros(Mx+2,NoTau);
  u(2:Mx+1,j) = B\f;
  
  if !isAsian
	u(inx-1,2:NoTau) = u(inx-1,NoTau+1);
  endif
  u(1,2:NoTau)=left_cond(m);
  lastu = u;

  alphap = -sigma(2:Mx,Mt+1-m).^2/2/dx^2 +(sigma(2:Mx,Mt+1-m).^2/2-(r0-r1))/2/dx;
  alpham = -sigma(2:Mx,Mt+1-m).^2/2/dx^2 -(sigma(2:Mx,Mt+1-m).^2/2-(r0-r1))/2/dx;
  beta = sigma(2:Mx,Mt+1-m).^2/dx^2+r0;
  A = spdiags([alpham, beta, alphap], -1:1, Mx, Mx);
  # compute matrices for the theta scheme
  theta = 0.5;
  B = speye(Mx,Mx) + theta*dt*A;
  C = speye(Mx,Mx) - (1-theta)*dt*A;
  
  if m == Mt-1
	 V2 = lastu(:,NoTau+1);
  endif
end

x = exp(x);
V1=lastu(:,NoTau+1);
endfunction

function [res] = CalculatePriceGreeksParisian(F_bid, F_ask, barrier,day_hat, strike,
			   issue_date,expire_date,PPO,OSO,price_type, \
                           barrier_type, payoff_type,isAsian)
  
global dsigma;


[x,V1,V2,F,S0,dx,dt,DF_d,DF_f] = Parisian_out(F_bid, F_ask, barrier,day_hat, strike,
			   issue_date,expire_date,PPO,OSO,price_type, \
                           barrier_type, payoff_type,isAsian,0);
  %% Price
  res(1) = interp1(x, V1, S0);
  
  %% Greeks

  %% spot delta
  res(2) = \ 
  ( interp1(x, V1, S0 + dx) \ 
	-interp1(x, V1, S0 - dx) ) / (2 * dx);
  
  %% forward delta
  res(3) = \ 
  ( interp1(x, V1, (F + dx) * DF_d / DF_f) \ 
	-interp1(x, V1, (F - dx) * DF_d / DF_f) ) / (2 * dx);
  
  %% spot gamma
  res(4) = \ 
  ( interp1(x, V1, S0 + dx) \ 
	-2 * interp1(x, V1, S0) \
	+interp1(x, V1, S0 - dx)) / dx^2;
  
  %% forward gamma
  res(5) = \ 
  ( interp1(x, V1, (F + dx) * DF_d / DF_f) \ 
	-2 * interp1(x, V1, S0) \
	+interp1(x, V1, (F - dx) * DF_d / DF_f) ) / dx^2;
  
  %% theta
  res(6) = (interp1(x, V2, S0) - res(1))/ dt ;

  %% vega
 
  [x2,V12,V22,F2,S02,dx2,dt2,DF_d2,DF_f2] = Parisian_out(F_bid, F_ask, barrier,day_hat, strike,
			   issue_date,expire_date,PPO,OSO,price_type, \
                           barrier_type, payoff_type,isAsian,1);
  res(7) = ( res(1) - interp1(x2, V12, S02) ) / (2 * dsigma);

endfunction
