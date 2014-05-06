1;
function [result_DM_doamcall] = DM_doamcall(F_bid, F_ask, barrier, strike, monitoring_dates, issue_date,expire_date,PPO,OSO,type)
 
 %% Market globals
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 
 %% PDE grid globals
 global Mx;
 global Mt;
 
 tau = year_frac(issue_date,expire_date,DCC);

 DF_d_bid = DF(issue_date,{expire_date},DSQ_Bid)(1);
 DF_d_ask = DF(issue_date,{expire_date},DSQ_Ask)(1);
 DF_f_bid = DF(issue_date,{expire_date},DSB_Bid)(1);
 DF_f_ask = DF(issue_date,{expire_date},DSB_Ask)(1);
 
 # data preparation for discounting with real dates of cash flows
 issue_mod = day_shift2(issue_date,FC_DOM, PPO);
 expire_mod = day_shift2(expire_date,FC_DOM, OSO);
 tau_mod = year_frac(issue_mod,expire_mod,DCC);
 
 # modified discount factors for real dates of cash flows
 DF_d_bid_mod = DF(issue_mod,{expire_mod},DSQ_Bid)(1);
 DF_d_ask_mod = DF(issue_mod,{expire_mod},DSQ_Ask)(1);
 
 if (type == "bid")
   DF_dG = DF_d_ask_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
  else
   DF_dG = DF_d_bid_mod;
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
 endif


 % monitoring dates adjustment
 [t, monitoring_ind] = adj_time_grid (tau, Mt, issue_date, expire_date, monitoring_dates);
  
 %% parameters of BS equation
 r0 = -log(DF_d)/tau;
 r1 = -log(DF_f)/tau;
 
 S0 = F * DF_d / DF_f;
 xmin = 0;
 xmax = set_xmax(S0, tau, vol(0, S0, strike, tau));
 x = linspace(xmin, xmax, Mx + 1);
 t = linspace(0, tau, Mt + 1);
 dx = (xmax - xmin)/Mx;
 dt = tau / Mt;
 sigma = vol(t, x, strike, tau);

 term_cond = max(x - strike, 0);
 left_cond = zeros(Mt + 1, 1);
 right_cond = (xmax - strike) * ones(Mt + 1, 1);
 
 left_bound = barrier * monitoring_ind;
 right_bound = xmax * ones(Mt + 1, 1);
 
 
 dsigma = 0.001;
 result_DM_doamcall = calc_price_greeks (r0, r1, sigma, term_cond, left_cond, right_cond, \
				   left_bound, right_bound, xmin, xmax, tau, Mx, Mt, \
				   S0, dsigma);
 
 

endfunction


function res = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, \
			   right_cond, left_bound, right_bound, xmin, xmax, T, \
			   Mx, Mt, res_t_size)
	 
	 if nargin < 14
		res_t_size = 1; %% result for t=0 only
	 endif
	 
	 res = -ones(res_t_size, Mx + 1);
	 
	 V = term_cond';
	 x = linspace(xmin, xmax, Mx + 1)';
	 dt = T / Mt;
	 dx = (xmax - xmin)/Mx;

	 for k = Mt:-1:1

	   ## Auxillary vectors for numerical scheme
	   is = xmin/dx + (0:Mx)';
	   r = r0 - r1;
	   A = .5*(sigma(:,k).^2.*is.^2-r*is)*dt;
	   B = -(sigma(:,k).^2.*is.^2 + r0)*dt;
	   C = .5*(sigma(:,k).^2.*is.^2 + r*is)*dt;

	   ## Only Crank - Nicholson so far
	   M1 = -spdiags([.5*C(2:Mx), .5*B(2:Mx)-1, .5*A(2:Mx)], [-1 0 1], Mx-1, Mx-1)';
	   M2 = -spdiags([-.5*C(2:Mx), -.5*B(2:Mx)-1, -.5*A(2:Mx)], [-1 0 1], Mx-1, Mx-1)';
	   d = zeros(Mx - 1, 1);

	   ## Applying boundary condition
	   d(1) = .5*A(1) * (left_cond(k) + left_cond(k+1));
	   d(end) = .5*C(end) * (right_cond(k) + right_cond(k+1));
	   V(2:Mx) = M1\(M2*V(2:Mx) + d);
	   V = V .* (x > left_bound(k) - dx / 2 & x < right_bound(k) + dx/2);
	
	   res(min(k, res_t_size),:) = V;
   
	endfor

endfunction


function result_vol = vol (t, S, K, T)
  result_vol = repmat(0.2, length(S), length(t));
#ImpVol(strike,issue_date,expire_date);
endfunction

function xmax = set_xmax(S0, T, sigma) 
  xmax = S0 * exp(6 * sigma * sqrt(T)); %% 6 sigmas
endfunction

%% Function adjusting time grid to barrier monitoring dates
function [t_grid, monitoring_ind] = adj_time_grid (T, Mt, issue_date, expire_date, monitoring_dates)
 
 global DCC;

 monitoring_t = zeros(size(monitoring_dates));

 for i = 1:length(monitoring_dates)
  monitoring_t(i) = year_frac(issue_date, monitoring_dates(i), DCC);
 endfor

 MF_precision = 100000;
 Mt_adj = round(MF_precision * T / gcd(round(MF_precision * monitoring_t)(2), round(MF_precision * monitoring_t)(3)));
 Mt = ceil(Mt / Mt_adj) * Mt_adj;
 dt = T / Mt;
 monitoring_k = round(monitoring_t / dt) + 1;
 monitoring_ind = zeros(Mt + 1, 1); 
 
 monitoring_ind(monitoring_k) = 1;
 t_grid = linspace(0, T, Mt + 1);
 
endfunction

function res = calc_price_greeks (r0, r1, sigma, term_cond, left_cond, right_cond, \
				   left_bound, right_bound, xmin, xmax, T, Mx, Mt, \
				   S0, dsigma)


  res = zeros(1, 7);
  dx = (xmin - xmax) / Mx;
  dt = T / Mt;
  DF_d = exp(-r0 * T);
  DF_f = exp(-r1 * T);
  F = S0 * DF_f / DF_d;

  %% Solving BS PDE
  x = linspace(xmin, xmax, Mx + 1);
  V = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
					left_bound, right_bound, xmin, xmax, T, Mx, Mt, 2);
  
  
  %% Price
  res(1) = interp1(x, V(1, :), S0);
  
  %% Greeks

  %% spot delta
  res(2) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-interp1(x, V(1, :), S0 - dx) ) / (2 * dx);
  
  %% forward delta
  res(3) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / (2 * dx);
  
  %% spot gamma
  res(4) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), S0 - dx)) / dx^2;
  
  %% forward gamma
  res(5) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / dx^2;
  
  %% theta
  res(6) = (interp1(x, V(2, :), S0) - res(1))/ dt ;

  %% vega
  V1 = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond, left_cond, right_cond, left_bound, right_bound, xmin, xmax, T, Mx, Mt);
  V2 = solve_BS_PDE (r0, r1, sigma - dsigma, term_cond, left_cond, right_cond, left_bound, right_bound, xmin, xmax, T, Mx, Mt);
  
  res(7) = \ 
  ( interp1(x, V1, S0) \ 
	-interp1(x, V2, S0) ) / (2 * dsigma);
  
		 
endfunction

