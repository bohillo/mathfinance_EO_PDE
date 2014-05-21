1;
wroblewski_ap_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%%% Functions for option valuation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = DM_out(F_bid, F_ask, barrier, strike,
			   monitoring_dates,
			   issue_date,expire_date,PPO,OSO,price_type, \
                           barrier_type, payoff_type)
  %% PDE grid globals
  global Mx;
  global Mt;

  if (price_type == "bid")
    F = F_bid;
  else
    F = F_ask;
  endif
  
  [tau, tau_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, price_type);
  % monitoring dates adjustment
  [t, monitoring_ind] = adj_time_grid (tau, Mt, issue_date, expire_date, monitoring_dates);
  Mt = length(t) - 1;
  %% parameters of BS equation 
  r0 = -log(DF_d)/tau;
  r1 = -log(DF_f)/tau;
  S0 = F * DF_d / DF_f;
  xmin = 0;
  xmax = set_xmax(S0, tau, vol(0, S0, strike, tau));
  x = linspace(xmin, xmax, Mx + 1);
  dx = (xmax - xmin)/Mx;
  dt = tau / Mt;
  sigma = vol(t, x, strike, tau);

  if strcmp(payoff_type,"put")
    left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
    right_cond = zeros(Mt + 1, 1);
    term_cond = max(strike - x, 0); 
  elseif strcmp(payoff_type,"call")
    term_cond = max(x - strike, 0);
    left_cond = zeros(Mt + 1, 1);
    right_cond = (xmax - strike) * ones(Mt + 1, 1);
  else
    left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
    right_cond = zeros(Mt + 1, 1);
    term_cond = max(x - strike, 0);   
  endif

  if strcmp(barrier_type,"up") 
    left_bound =  zeros(Mt + 1, 1);
    right_bound = barrier * monitoring_ind + xmax * !monitoring_ind;
  elseif strcmp(barrier_type, "down")
    right_bound =  xmax * ones(Mt + 1, 1);
    left_bound = barrier * monitoring_ind;
  else
    left_bound =  zeros(Mt + 1, 1);
    right_bound = xmax * ones(Mt + 1, 1);
  endif
 
 
  result = calc_price_greeks (r0, r1, sigma, term_cond, left_cond, right_cond, \
			      left_bound, right_bound, xmin, xmax, tau, Mx, Mt, \
			      S0);

  result(1) = result(1) * DF_dG / DF_d;
endfunction

function [result] = DM_in(F_bid, F_ask, barrier, strike,
				      monitoring_dates,
				      issue_date,expire_date,PPO,OSO,price_type, \
                                      barrier_type, payoff_type)

 if strcmp(payoff_type, "call")
   result = \
   call(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,price_type) - \
   DM_out(F_bid, F_ask, barrier, strike,
				      monitoring_dates,
				      issue_date,expire_date,PPO,OSO,price_type, \
                                      barrier_type, payoff_type);
 elseif strcmp(payoff_type, "put") 
   result = \
   put(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,price_type) - \
   DM_out(F_bid, F_ask, barrier, strike,
				      monitoring_dates,
				      issue_date,expire_date,PPO,OSO,price_type, \
                                      barrier_type, payoff_type);
 else
	error("Invalid payoff type");
 endif
endfunction


function [result] = DoubleKO(F_bid, F_ask, Lbarrier, Ubarrier, strike,
			   monitoring_dates,
			   issue_date,expire_date,PPO,OSO,price_type, \
			   payoff_type)

  %% PDE grid globals
  global Mx;
  global Mt;

  if (price_type == "bid")
    F = F_bid;
  else
    F = F_ask;
  endif
  
  [tau, tau_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, price_type);
  % monitoring dates adjustment
  [t, monitoring_ind] = adj_time_grid (tau, Mt, issue_date, expire_date, monitoring_dates);
  Mt = length(t) - 1;
  %% parameters of BS equation 
  r0 = -log(DF_d)/tau;
  r1 = -log(DF_f)/tau;
  S0 = F * DF_d / DF_f;

  if isempty(monitoring_dates)
    xmin = Lbarrier;
    xmax = Ubarrier;
  else
    xmin = 0;
    xmax = set_xmax(S0, tau, vol(0, S0, strike, tau));
  endif

  x = linspace(xmin, xmax, Mx + 1);
  dx = (xmax - xmin)/Mx;
  dt = tau / Mt;
  sigma = vol(t, x, strike, tau);

  if strcmp(payoff_type,"put")
    left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
    right_cond = zeros(Mt + 1, 1);
    term_cond = max(strike - x, 0); 
  elseif strcmp(payoff_type,"call")
    term_cond = max(x - strike, 0);
    left_cond = zeros(Mt + 1, 1);
    right_cond = (xmax - strike) * ones(Mt + 1, 1);
  else
    left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
    right_cond = zeros(Mt + 1, 1);
    term_cond = max(x - strike, 0);   
  endif

  if isempty(monitoring_dates)
    right_cond = zeros(Mt + 1, 1);
    left_cond = zeros(Mt + 1, 1);
  endif

  left_bound = Lbarrier * monitoring_ind + xmin * !monitoring_ind;
  right_bound = Ubarrier * monitoring_ind + xmax * !monitoring_ind;
 
  result = calc_price_greeks (r0, r1, sigma, term_cond, left_cond, right_cond, \
			      left_bound, right_bound, xmin, xmax, tau, Mx, Mt, \
			      S0);
  
  result(1) = result(1) * DF_dG / DF_d;
endfunction

function [result] = KIKO(F_bid, F_ask, Lbarrier, Ubarrier, strike,
			 monitoring_dates,
			 issue_date,expire_date,PPO,OSO,price_type, \
			 payoff_type)

  dbl_result = DoubleKO(F_bid, F_ask, Lbarrier, Ubarrier, strike,
			monitoring_dates,
			issue_date,expire_date,PPO,OSO,price_type, \
			payoff_type);

  if !isempty(monitoring_dates)
    result = dbl_result - DM_out(F_bid, F_ask, Ubarrier, strike,
			monitoring_dates,
			issue_date,expire_date,PPO,OSO,price_type, \
			"up", payoff_type);
  else
    result = dbl_result - uoamcall(F_bid,F_ask,Ubarrier,strike,issue_date,expire_date,PPO,OSO,price_type);
  endif 

endfunction


function [result] = Window_out(F_bid, F_ask, barrier, strike, \
			       issue_date, window_start_date, window_end_date, \
			       expire_date,PPO,OSO,price_type, barrier_type, payoff_type)
  %% PDE grid globals
  global Mx;
  global Mt;

  %% For Vega computation
  global dsigma;

  if (price_type == "bid")
    F = F_bid;
  else
    F = F_ask;
  endif
  
  Mt_tot = Mt;
  
  [tau1, tau_mod1, DF_dG1, DF_d1, DF_f1] = prepare_data(window_end_date, expire_date, 0, PPO, price_type);
  [tau2, tau_mod2, DF_dG2, DF_d2, DF_f2] = prepare_data(window_start_date, window_end_date, 0, 0, price_type);
  [tau3, tau_mod3, DF_dG3, DF_d3, DF_f3] = prepare_data(issue_date, window_start_date, OSO, 0, price_type);
  [tau, tau_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, price_type);

  S0 = F * DF_d / DF_f;
  
  %% parameters of BS equation 
  % 1st step

  xmin1 = 0;
  xmax1 = set_xmax(S0, tau, vol(0, S0, strike, tau));
  x1 = linspace(xmin1, xmax1, Mx + 1);
  x= x1;
  dx = (xmax1 - xmin1)/Mx;
  
  if strcmp(payoff_type,"put")
    term_cond = max(strike - x, 0); 
  elseif strcmp(payoff_type,"call")
    term_cond = max(x - strike, 0);
  else
	error("Invalid payoff type");
  endif
  
  term_cond_vega = term_cond;

  if (tau1 > 0)
	Mt = ceil(Mt_tot * tau1 / tau);
	t = linspace(0, tau1, Mt + 1);
	r0 = -log(DF_d1)/tau1;
	r1 = -log(DF_f1)/tau1;
	dt = tau1 / Mt;
	sigma = vol(t, x, strike, tau1);

	if strcmp(payoff_type,"put")
      left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
      right_cond = zeros(Mt + 1, 1);
  	elseif strcmp(payoff_type,"call")
      left_cond = zeros(Mt + 1, 1);
      right_cond = (xmax1 - strike) * ones(Mt + 1, 1);
	else
	  error("Invalid payoff type");
	endif

	left_bound =  zeros(Mt + 1, 1);
	right_bound = xmax1 * ones(Mt + 1, 1);
	
	term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, \
							  tau1, Mx, Mt);
	term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau1, Mx, Mt);
  endif

  % 2nd step
  Mt = ceil(Mt_tot * tau2 / tau);
  t = linspace(0, tau2, Mt + 1);
  r0 = -log(DF_d2)/tau2;
  r1 = -log(DF_f2)/tau2;
  xmin = xmin1;
  xmax = xmax1;

  if strcmp(barrier_type,"up") 
	xmax = barrier;
  elseif strcmp(barrier_type, "down")
	xmin = barrier;
  else
    error("Invalid barrier type");
  endif
  
  x = x(x1 >= xmin  & x1 <= xmax);
  dt = tau2 / Mt;
  sigma = vol(t, x, strike, tau2);

  if strcmp(payoff_type,"put")
    left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
    right_cond = zeros(Mt + 1, 1);
  elseif strcmp(payoff_type,"call")
    left_cond = zeros(Mt + 1, 1);
    right_cond = (xmax - strike) * ones(Mt + 1, 1);
  else
	error("Invalid payoff type");
  endif

  if strcmp(barrier_type,"up") 
    right_cond = zeros(Mt + 1, 1);
    left_bound =  zeros(Mt + 1, 1);
    right_bound = barrier * ones(Mt + 1, 1);
  elseif strcmp(barrier_type, "down")
    left_cond = zeros(Mt + 1, 1);
    right_bound =  xmax * ones(Mt + 1, 1);
    left_bound = barrier * ones(Mt + 1, 1);
  else
    error("Invalid barrier type");
  endif

  term_cond = term_cond(x1>= xmin  & x1 <= xmax);
  term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
					left_bound, right_bound, xmin, xmax, tau2, \
					length(term_cond) - 1, Mt, 2);

  term_cond_vega = term_cond_vega(x1>= xmin  & x1 <= xmax);
  term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
					left_bound, right_bound, xmin, xmax, tau2, \
					length(term_cond) - 1, Mt);

  % 3rd step
  if (tau3 > 0)
	Mt = ceil(Mt_tot * tau3 / tau);
	t = linspace(0, tau3, Mt + 1);
	r0 = -log(DF_d3)/tau3;
	r1 = -log(DF_f3)/tau3;
	x = linspace(xmin1, xmax1, Mx + 1);

	dt = tau3 / Mt;
	sigma = vol(t, x, strike, tau3);

	if strcmp(payoff_type,"put")
      left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
      right_cond = zeros(Mt + 1, 1);
	elseif strcmp(payoff_type,"call")
      left_cond = zeros(Mt + 1, 1);
      right_cond = (xmax - strike) * ones(Mt + 1, 1);
	else
	  error("Invalid payoff type");
	endif

	left_bound =  zeros(Mt + 1, 1);
	right_bound = xmax1 * ones(Mt + 1, 1);

	term_cond1 = zeros(1, Mx + 1); 
	term_cond1(x >= xmin & x <= xmax) = term_cond(1,:);
	term_cond = term_cond1;
	term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau3, Mx, Mt);
  
	term_cond1 = zeros(1, Mx + 1); 
	term_cond1(x >= xmin & x <= xmax) = term_cond_vega;
	term_cond_vega = term_cond1;
	term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau3, Mx, Mt);
  endif
  
  V = term_cond;
  V_vega = term_cond_vega;
  
  %% Price
  result(1) = interp1(x, V(1, :), S0);
  
  %% Greeks

  %% spot delta
  result(2) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-interp1(x, V(1, :), S0 - dx) ) / (2 * dx);
  
  %% forward delta
  result(3) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / (2 * dx);
  
  %% spot gamma
  result(4) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), S0 - dx)) / dx^2;
  
  %% forward gamma
  result(5) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / dx^2;
  
  %% theta
  result(6) = (interp1(x, V(2, :), S0) - result(1))/ dt ;
  
  result(7) = \ 
  ( interp1(x, V_vega, S0) \ 
	-interp1(x, V(1,:), S0) ) /  dsigma;
  
  result(1) = result(1) * DF_dG / DF_d;
endfunction

function [result] = Window_in(F_bid, F_ask, barrier, strike, \
			   issue_date, window_start_date, window_end_date, \
			   expire_date,PPO,OSO,price_type, barrier_type, payoff_type)

 if strcmp(payoff_type, "call")
   result = \
   call(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,price_type) - \
   Window_out(F_bid, F_ask, barrier, strike, issue_date, window_start_date, window_end_date, \
			   expire_date,PPO,OSO,price_type, barrier_type, payoff_type);
 elseif strcmp(payoff_type, "put") 
   result = \
   put(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,price_type) - \
   Window_out(F_bid, F_ask, barrier, strike, issue_date, window_start_date, window_end_date, \
			   expire_date,PPO,OSO,price_type, barrier_type, payoff_type);
 else
	error("Invalid payoff type");
 endif
		 
endfunction 

function [result] = Window_DoubleKO(F_bid, F_ask, Lbarrier, Ubarrier, strike,\
				    issue_date, window_start_date, window_end_date, \
				    expire_date,PPO,OSO,price_type, payoff_type)
  %% PDE grid globals
  global Mx;
  global Mt;

  %% For Vega computation
  global dsigma;

  if (price_type == "bid")
    F = F_bid;
  else
    F = F_ask;
  endif
  
  Mt_tot = Mt;
  
  [tau1, tau_mod1, DF_dG1, DF_d1, DF_f1] = prepare_data(window_end_date, expire_date, 0, PPO, price_type);
  [tau2, tau_mod2, DF_dG2, DF_d2, DF_f2] = prepare_data(window_start_date, window_end_date, 0, 0, price_type);
  [tau3, tau_mod3, DF_dG3, DF_d3, DF_f3] = prepare_data(issue_date, window_start_date, OSO, 0, price_type);
  [tau, tau_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, price_type);

  S0 = F * DF_d / DF_f;
  
  %% parameters of BS equation 
  % 1st step

  xmin1 = 0;
  xmax1 = set_xmax(S0, tau, vol(0, S0, strike, tau));
  x1 = linspace(xmin1, xmax1, Mx + 1);
  x= x1;
  dx = (xmax1 - xmin1)/Mx;
  
  if strcmp(payoff_type,"put")
    term_cond = max(strike - x, 0); 
  elseif strcmp(payoff_type,"call")
    term_cond = max(x - strike, 0);
  else
	error("Invalid payoff type");
  endif
  
  term_cond_vega = term_cond;

  if (tau1 > 0)
	Mt = ceil(Mt_tot * tau1 / tau);
	t = linspace(0, tau1, Mt + 1);
	r0 = -log(DF_d1)/tau1;
	r1 = -log(DF_f1)/tau1;
	dt = tau1 / Mt;
	sigma = vol(t, x, strike, tau1);

	if strcmp(payoff_type,"put")
      left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
      right_cond = zeros(Mt + 1, 1);
  	elseif strcmp(payoff_type,"call")
      left_cond = zeros(Mt + 1, 1);
      right_cond = (xmax1 - strike) * ones(Mt + 1, 1);
	else
	  error("Invalid payoff type");
	endif

	left_bound =  zeros(Mt + 1, 1);
	right_bound = xmax1 * ones(Mt + 1, 1);
	
	term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, \
							  tau1, Mx, Mt);
	term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau1, Mx, Mt);
  endif

  % 2nd step
  Mt = ceil(Mt_tot * tau2 / tau);
  t = linspace(0, tau2, Mt + 1);
  r0 = -log(DF_d2)/tau2;
  r1 = -log(DF_f2)/tau2;
  xmin = Lbarrier;
  xmax = Ubarrier;

  x = x(x1 >= xmin  & x1 <= xmax);
  dt = tau2 / Mt;
  sigma = vol(t, x, strike, tau2);


  left_cond = zeros(Mt + 1, 1);
  right_cond = zeros(Mt + 1, 1);
  left_bound =  Lbarrier * ones(Mt + 1, 1);
  right_bound = Ubarrier * ones(Mt + 1, 1);
  
  term_cond = term_cond(x1>= xmin  & x1 <= xmax);
  term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
					left_bound, right_bound, xmin, xmax, tau2, \
					length(term_cond) - 1, Mt, 2);

  term_cond_vega = term_cond_vega(x1>= xmin  & x1 <= xmax);
  term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
					left_bound, right_bound, xmin, xmax, tau2, \
					length(term_cond) - 1, Mt);

  % 3rd step
  if (tau3 > 0)
	Mt = ceil(Mt_tot * tau3 / tau);
	t = linspace(0, tau3, Mt + 1);
	r0 = -log(DF_d3)/tau3;
	r1 = -log(DF_f3)/tau3;
	x = linspace(xmin1, xmax1, Mx + 1);

	dt = tau3 / Mt;
	sigma = vol(t, x, strike, tau3);

	if strcmp(payoff_type,"put")
      left_cond = strike * ones(Mt + 1, 1) .* exp(-r0 * (tau - t'));
      right_cond = zeros(Mt + 1, 1);
	elseif strcmp(payoff_type,"call")
      left_cond = zeros(Mt + 1, 1);
      right_cond = (xmax - strike) * ones(Mt + 1, 1);
	else
	  error("Invalid payoff type");
	endif

	left_bound =  zeros(Mt + 1, 1);
	right_bound = xmax1 * ones(Mt + 1, 1);

	term_cond1 = zeros(1, Mx + 1); 
	term_cond1(x >= xmin & x <= xmax) = term_cond(1,:);
	term_cond = term_cond1;
	term_cond = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau3, Mx, Mt);
  
	term_cond1 = zeros(1, Mx + 1); 
	term_cond1(x >= xmin & x <= xmax) = term_cond_vega;
	term_cond_vega = term_cond1;
	term_cond_vega = solve_BS_PDE (r0, r1, sigma + dsigma, term_cond_vega, left_cond, right_cond, \
							  left_bound, right_bound, xmin1, xmax1, tau3, Mx, Mt);
  endif
  
  V = term_cond;
  V_vega = term_cond_vega;
  
  %% Price
  result(1) = interp1(x, V(1, :), S0);
  
  %% Greeks

  %% spot delta
  result(2) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-interp1(x, V(1, :), S0 - dx) ) / (2 * dx);
  
  %% forward delta
  result(3) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / (2 * dx);
  
  %% spot gamma
  result(4) = \ 
  ( interp1(x, V(1, :), S0 + dx) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), S0 - dx)) / dx^2;
  
  %% forward gamma
  result(5) = \ 
  ( interp1(x, V(1, :), (F + dx) * DF_d / DF_f) \ 
	-2 * interp1(x, V(1, :), S0) \
	+interp1(x, V(1, :), (F - dx) * DF_d / DF_f) ) / dx^2;
  
  %% theta
  result(6) = (interp1(x, V(2, :), S0) - result(1))/ dt ;
  
  result(7) = \ 
  ( interp1(x, V_vega, S0) \ 
	-interp1(x, V(1,:), S0) ) /  dsigma;
  
  result(1) = result(1) * DF_dG / DF_d;
endfunction

function [result] = Window_KIKO(F_bid, F_ask, Lbarrier, Ubarrier, strike, \ 
			 issue_date, window_start_date, window_end_date, expire_date,PPO,OSO,price_type, \
			 payoff_type)

  dbl_result = Window_DoubleKO(F_bid, F_ask, Lbarrier, Ubarrier, strike, \
			issue_date, window_start_date, window_end_date, expire_date,PPO,OSO,price_type, \
			payoff_type);

    result = dbl_result - Window_out(F_bid, F_ask, Ubarrier, strike, \
			issue_date, window_start_date, window_end_date, expire_date,PPO,OSO,price_type, \
			"up", payoff_type);

endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%%% Auxilary functions  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = ImpVol(strike,issue_date,expire_date)
  global DCC;
  result = vol(0, strike, strike, year_frac(issue_date, expire_date, DCC) );
endfunction

                     
function [tau, tau_mod, DF_dG, DF_d, DF_f] = prepare_data(issue_date, expire_date, OSO, PPO, type)  
   %% Market globals
   global DCC;
   global FC_DOM;
   global DSQ_Bid;
   global DSQ_Ask;
   global DSB_Bid;
   global DSB_Ask;

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
	 DF_d = DF_d_bid;
	 DF_f = DF_f_ask;
   else
	 DF_dG = DF_d_bid_mod;
	 DF_d = DF_d_ask;
	 DF_f = DF_f_bid;
   endif
 
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

	   ## Crank - Nicholson
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
endfunction

function xmax = set_xmax(S0, T, sigma) 
  xmax = S0 * exp(6 * sigma * sqrt(T)); %% 6 sigmas
endfunction

%% Function adjusting time grid to barrier monitoring dates
function [t_grid, monitoring_ind] = adj_time_grid (T, Mt, issue_date, expire_date, monitoring_dates)
 
 global DCC;

 if isempty(monitoring_dates) 
   t_grid = linspace(0, T, Mt + 1);
   monitoring_ind = ones(1, Mt + 1);
   return;
 endif

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
				   S0)

  global dsigma;
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

  res(7) = \ 
  ( interp1(x, V1, S0) \ 
	-res(1) ) / dsigma;
  
		 
endfunction

