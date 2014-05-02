% TO JEST WERSJA DO UZYTKU POD OCTAVE
% NA KONCU SA FUNKCJE KTORE ZAPISUJA DO PLIKOW
% INFORMACJE POTRZEBNE DO URUCHOMIENIA CZESCI "PRICES"


% powierzchnia zmiennosci implikowanej Joanna Pieniak, Piotr Wiazecki 


% Standardowe oznaczenia:
% S_0 - wartosæ poczatkowa 
% K - strike
% f = FX_rate * DF_f/DF_d;
% tau - u³amek roku
% DF_f -   czynnik dyskontowy waluty bazowej
% DF_d -   czynnik dyskontowy waluty niebazowej
% poprawki 07.09.13
% poprawiono bledne wzory na kurs forward
%

1;
 
############################ GENEROWANIE fx_vol.m ################################
% funkcja 1 --------------------------------------------------------------------
% generowanie fx_vol.m

function  market_vol_gen()
global VOL_DCC; 
global VOL_BDA; 
global VOL_EMA; 
global FX_rate; 
global vol_fx_rate; 
global vol_start_date;
global FX_VOL;
global premium_curr;
global FC_DOM;
global FC_FOR;
global BAS_CURR;
global QUO_CURR;
global ST_BAS_QUO; 
 
   
  premium_curr =  PREM_CURR(BAS_CURR, QUO_CURR);
   
  fid = fopen("market_vol.m", "w");
  fprintf(fid, "premium_curr = \"%s\";\n\n",  premium_curr);
  fprintf(fid, "vol_start_date = \"%s\";\n", vol_start_date);
  fprintf(fid, "vol_fx_rate = %f;\n",  vol_fx_rate);
  fprintf(fid, "FX_VOL_REC = {");
  
  if (strcmp(ST_BAS_QUO, "true") ==1)
     FC_VOL = FC_DOM;
  else 
     FC_VOL = FC_FOR;
  endif
         
  
  no_rows_FX_VOL = rows(FX_VOL);  
  for i = 1 : rows(FX_VOL)
    fprintf(fid, "\"%s\",",  FX_VOL{i, 1});
    [pay_dt, pay_fr] = cpn_dt(vol_start_date, FX_VOL{i, 1}, "DEP", FC_VOL, VOL_DCC, VOL_BDA, 0, VOL_EMA, 0); 
    fprintf(fid, "\"%s\",", pay_dt{2});
    fprintf(fid, "%f,", FX_VOL{i, 2});
    fprintf(fid, "%f,", FX_VOL{i, 3});
    fprintf(fid, "%f,", FX_VOL{i, 4});
	 
    if index(FX_VOL{i, 1}, "Y") > 0
      fprintf(fid, "\"atm_fwd\",");
    else
      fprintf(fid, "\"delta_neutral\",");
    end
	
    days = day_diff(vol_start_date, pay_dt{2}, "ACT");
    delta_conv = delta_convention(BAS_CURR, QUO_CURR, days);
    fprintf(fid, "\"%s\"", delta_conv);
    if i < no_rows_FX_VOL
    	fprintf(fid, ";", delta_conv);
    endif
  end
  
  
  fprintf(fid, "};\n");
  fclose(fid);
  
endfunction

% funkcja 2 ------------------------------------------------------------------------------
% funkcja PREM_CURR oblicza walute premium
% input: CURR_FOR - waluta obca, CURR_DOM - waluta krajowa
% output: waluta premium

function result =  PREM_CURR(CURR_FOR, CURR_DOM)

	order = {"USD", "EUR", "GBP", "AUD", "NZD", "CAD", "CHF", "NOK", "SEK", "DKK", "CZK", "PLN", "TRY", "MXN", "JPY"};

	for i = 1 : columns(order)
		if strcmp(order(i), CURR_FOR) == 1
			wsk_c = i;
		end
		if strcmp(order(i), CURR_DOM) == 1
			wsk_d = i;
		end
	end

	if wsk_c > wsk_d
		PREM_CURR = CURR_DOM;
		else PREM_CURR = CURR_FOR;
	end 

	result = PREM_CURR;

endfunction

% funkcja 3  --------------------------------------------------------------------
% funkcja oblicza konwencje delty dla danej pary walutowej
% input: 
%   CURR_FOR - waluta obca,
%   CURR_DOM - waluta krajowa
%   tau - czas do zapadalnosci
% output: konwencja delty

function result = delta_convention(BAS_CURR, QUO_CURR, tau)
global premium_curr;


% spot czy fwd
	OECD = {"USD", "EUR", "JPY", "GBP", "AUD", "NZD", "CAD", "CHF", "NOK", "SEK", "DKK"};
	for i = 1 : columns(OECD)
		a(i) = strcmp(OECD(i), QUO_CURR);
		b(i) = strcmp(OECD(i), BAS_CURR);
	end

	w = a + b;
	if (sum(w) == 2 && tau < 366) 
		delta_conv(1) = 0;
		else delta_conv(1) = 1;
	end

% premium adjusted czy unadjusted 

	wsk_c = 0;
	if strcmp(premium_curr, QUO_CURR) == 1
			wsk_c = 1;
	end
	if wsk_c == 1
		delta_conv(2) = 0;
	else 
		delta_conv(2) = 1;
	end
	% delta = (0,0) - spot, unadjusted
	% delta = (0,1) - spot, premium adj
	% delta = (1,0) - fwd, unadjusted
	% delta = (1,1) - fwd, premium adj

	conversions = {"SPOT_UN", "SPOT_PA"; "FORWARD_UN", "FORWARD_PA"};
	result = conversions{delta_conv(1)+1, delta_conv(2)+1};
endfunction

% function 4 ---------------------------------------------------------------
% opis create_vol_table(): wczytuje dane z pliku market_vol
%  przygotowuje dane (ceny wykonania, zmiennosci i wspó³czynniki dyskontowe) do pózniejszego u¿ycia


function create_vol_table()
  global smile_interp;
  global m_points;
  global volatility_grid;
  global K_grid;
  global T_grid;
  global fx_vol_data;
  
	source market_vol.m;
    vol_rec = FX_VOL_REC;
    vol_grid = zeros(rows(vol_rec),m_points);
    K = zeros(rows(vol_rec),m_points);
    fx_vol_data = zeros(rows(vol_rec),8);
    T_grid = zeros(1, rows(vol_rec));	 
	
	 for i = 1 : rows(vol_rec)
	   
	     sigma_BF2vol = BF2vol(vol_rec{i, 2}, vol_rec{i, 3}, vol_rec{i, 4}, vol_rec{i, 5}, vol_rec{i, 6}, vol_rec{i, 7});
		 
	     [strikes, sigmas, DFs, tau, limits] = VOL(vol_rec{i, 2}, vol_rec{i, 3}, vol_rec{i, 4}, sigma_BF2vol, vol_rec{i, 6}, vol_rec{i, 7});
         fx_vol_data(i,:) = [strikes(1), strikes(2), strikes(3), sigmas(1), sigmas(2), sigmas(3), DFs(1), DFs(2)];
         T_grid(i) = tau;
         K(i,:) = linspace(limits(1),limits(2),m_points);   
         for j = 1 : m_points    
             vol_grid(i,j) = strike_interpolate(strikes, sigmas, DFs, tau,  K(i,j), smile_interp);
         endfor  
    endfor
   
   
   volatility_grid = vol_grid;
   K_grid = K;
   
   fid = fopen("vol_tables.m", "w");
   
   fprintf(fid, "volatility_grid = [");
    
   for i = 1 : rows(vol_rec)
     for j = 1 : (m_points-1) 
      fprintf(fid, "%f,", volatility_grid(i,j));
     endfor
     fprintf(fid, "%f;\n", volatility_grid(i,m_points));
   endfor
   fprintf(fid, "];\n\n");
   
   fprintf(fid, "K_grid = [\n");
    
   for i = 1 : rows(vol_rec)
     for j = 1 : (m_points-1) 
      fprintf(fid, "%f,", K_grid(i,j));
     endfor
     fprintf(fid, "%f;\n", K_grid(i,m_points));
   endfor
   fprintf(fid, "];\n\n");
   
   fprintf(fid, "T_grid = [\n");
    
   for i = 1 : rows(vol_rec)
      fprintf(fid, "%f,", T_grid(i));
   endfor
   fprintf(fid, "];\n\n");
   
   fprintf(fid, "fx_vol_data = [\n");
    
   for i = 1 : rows(vol_rec)
     for j = 1 : 7 
      fprintf(fid, "%f,", fx_vol_data(i,j));
     endfor
     fprintf(fid, "%f;\n", fx_vol_data(i,8));
   endfor
   fprintf(fid, "];\n\n");
   
   fclose(fid);
  
endfunction

% funkcja 5 ------------------------------------------------
% opis VOL: zwraca wartosæ cen wykonania i zmiennosci dla deltaC, ATM, deltaP
% input: expiry, sigma_atm, sigma_rr, sigma_str - odczytane z rynku, ATM_conv, conv = {SPOT_UN, SPOT_PA, FORWARD_UN, FORWARD_PA}, premium_curr 
% output: strikes, sigmas, DFs, tau

function [strikes, sigmas, DFs, tau, limits] = VOL(expiry, sigma_atm, sigma_rr, sigma_str, ATM_conv, conv)

	global DF_QUOT;
   global DF_BASE; 
	global vol_start_date;
	global VOL_DCC;
	global vol_fx_rate;
    
	DF_f = DF(vol_start_date, {expiry}, DF_BASE); 
	DF_d = DF(vol_start_date, {expiry}, DF_QUOT);
	f = vol_fx_rate * DF_f/DF_d;
	S_0 = vol_fx_rate;
	tau = year_frac(vol_start_date, {expiry}, VOL_DCC);
   
	 
	sigma_atm /= 100;
	sigma_str /= 100;
	sigma_rr /= 100;

	sigma_C = sigma_atm + sigma_str + 0.5 * sigma_rr;
	sigma_P = sigma_atm + sigma_str - 0.5 * sigma_rr;

	if strcmp(conv, "SPOT_UN") == 1
		convention = spot_un(sigma_C, sigma_P, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv);
	elseif strcmp(conv, "SPOT_PA") == 1
		convention =  spot_pa(sigma_C, sigma_P, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv);
	elseif strcmp(conv, "FORWARD_UN") == 1
		convention = fwd_un(sigma_C, sigma_P, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv);
	elseif strcmp(conv, "FORWARD_PA") == 1 
		convention = fwd_pa(sigma_C, sigma_P, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv);
	end
   
	K_atm = convention(1);
	K_deltaC = convention(2);
	K_deltaP = convention(3);
	K_LC = convention(4);
	K_LP = convention(5);
	
	strikes = [K_deltaP, K_atm, K_deltaC];
	sigmas = [sigma_P, sigma_atm, sigma_C];
	DFs = [DF_f, DF_d];
	limits = [K_LP, K_LC];

endfunction

% funkcja 6 ------------------------------------------------
% opis BF2vol: wylicza poprawn¹ wartosæ sigma_BF2vol na podstawie danych rynkowych sigma_atm, sigma_rr, sigma_BF1vol, znajduj¹c miejsce zerowe
% funkcji strangle_value_difference opisanej poni¿ej
%
% input: expiry, sigma_atm, sigma_rr, sigma_BF1vol, ATM_conv, conv
% output: wartosæ sigma_BF2vol

function vol = BF2vol(expiry, sigma_atm, sigma_rr, sigma_BF1vol, ATM_conv, conv) 

    vol = fzero(@(sigma)strangle_value_difference(sigma, expiry, sigma_atm, sigma_rr, sigma_BF1vol, ATM_conv, conv),sigma_BF1vol); 
	
endfunction

% funkcja 7 ------------------------------------------------
% opis strangle_value_difference: oblicza ró¿nicê w wartosci broker's strangle wynikaj¹cej z podstawienia podanej na wejsciu wartosci sigma
% w miejsce sigma_BF2vol i skonstruowania krzywej volatility przy jej u¿yciu, oraz rynkowej wartosci sigma_BF1vol
% input: sigma, expiry, sigma_atm, sigma_rr, sigma_BF1vol, ATM_conv, conv

function diff = strangle_value_difference(sigma, expiry, sigma_atm, sigma_rr, sigma_BF1vol, ATM_conv, conv)

%   global DF_QUOT;
%   global DF_BASE; 
%	global vol_start_date;
%	global VOL_DCC;
	global vol_fx_rate;
	global RR_delta_val;
	global smile_interp;
	
	[strikes, sigmas, DFs, tau, limits] = VOL(expiry, sigma_atm, sigma_rr, sigma, ATM_conv, conv);
	DF_f = DFs(1);
	DF_d = DFs(2);
		
	f = vol_fx_rate * DF_f/DF_d;
	S_0 = vol_fx_rate;
	
	
	sigma_atm /= 100;
	sigma_rr /= 100;
	sigma_BF1vol /= 100;
	sigma /= 100;
   sigmaSTG1vol = sigma_BF1vol + sigma_atm;

	if strcmp(conv, "SPOT_UN") == 1
	
       K_deltaC_star = f * exp(-norminv( RR_delta_val/DF_f ) * sigmaSTG1vol * sqrt(tau) + 0.5 * sigmaSTG1vol^2 * tau);
	    K_deltaP_star = f * exp(norminv( RR_delta_val/DF_f ) * sigmaSTG1vol * sqrt(tau) + 0.5 * sigmaSTG1vol^2 * tau);
		
	elseif strcmp(conv, "SPOT_PA") == 1
	
	    e = 0.0001;
	    E = 100 * S_0;
	
	    K_deltaP_star = fzero(@(K)spot_pa_put(K, f, sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [e E]);
	    Kmin = fmax(@(K)spot_pa_call(K, f, sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [e E]);
	    K_deltaC_star = fzero(@(K)spot_pa_call(K, f,sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [Kmin K_deltaP_star]);
	
	elseif strcmp(conv, "FORWARD_UN") == 1
	
    	K_deltaC_star = f * exp(-norminv(RR_delta_val) * sigmaSTG1vol * sqrt(tau) + 0.5 * sigmaSTG1vol^2 * tau);
	    K_deltaP_star = f * exp(norminv(RR_delta_val) * sigmaSTG1vol * sqrt(tau) + 0.5 * sigmaSTG1vol^2 * tau);
		
	elseif strcmp(conv, "FORWARD_PA") == 1 
		
		e = 0.0001;
	    E = 100 * S_0;

	    K_deltaP_star = fzero(@(K)fwd_pa_put(K, f, sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [e E]);
	    Kmin = fmax(@(K)fwd_pa_call(K, f, sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [e E]);
	    K_deltaC_star = fzero(@(K)fwd_pa_call(K, f, sigmaSTG1vol, tau, DF_f, DF_d, RR_delta_val), [Kmin K_deltaP_star]);
		
	end

	sigma_call = strike_interpolate(strikes, sigmas, DFs, tau, K_deltaC_star, smile_interp);
	sigma_put = strike_interpolate(strikes, sigmas, DFs, tau, K_deltaP_star, smile_interp);
	
   Call_STG = C(f, K_deltaC_star, sigmaSTG1vol, tau, DF_d);
	Call_volsmile = C(f, K_deltaC_star, sigma_call, tau, DF_d);
	Put_STG = C(f, K_deltaP_star, sigmaSTG1vol, tau, DF_d) - S_0*DF_f + K_deltaP_star*DF_d;
	Put_volsmile = C(f, K_deltaP_star, sigma_put, tau, DF_d) - S_0*DF_f + K_deltaP_star*DF_d;
	
	diff = Call_STG + Put_STG - Call_volsmile - Put_volsmile;
endfunction

% funkcje 8-11 ------------------------------------------------
% opis spot_un, spot_pa, fwd_un, fwd_pa: obliczaja ceny wykonania K_atm, K_deltaC, K_deltaP 
% dla danych sigma_atm, sigma_deltaC, sigma_deltaP, dla poszczególnych konwencji 
% (w przypadku spot_pa i fwd_pa obliczane jest miejsce zerowe wyra¿eñ zdefiniwanych w pomocniczych funkcjach 
% fwd_pa_call, fwd_pa_put, spot_pa_call, spot_pa_put - opis poni¿ej, 
% korzystajac z wbudowanych w Octave funkcji fzero i napisanej funkcji fmax - opis poni¿ej)
%
% input: sigma_deltaC, sigma_deltaP, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv
% output: [K_atm, K_deltaC, K_deltaP]

function result = spot_un(sigma_deltaC, sigma_deltaP, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv)
   
   global RR_delta_val;
   global limit_delta;
	
		
	K_deltaC = f * exp(-norminv( RR_delta_val/DF_f ) * sigma_deltaC * sqrt(tau) + 0.5 * sigma_deltaC^2 * tau);
	K_deltaP = f * exp(norminv( RR_delta_val/DF_f ) * sigma_deltaP * sqrt(tau) + 0.5 * sigma_deltaP^2 * tau);
	K_atm = ATM_strike(S_0, f, sigma_atm, tau, ATM_conv, 1);
	K_LC = f * exp(-norminv( limit_delta/DF_f ) * sigma_deltaC * sqrt(tau) + 0.5 * sigma_deltaC^2 * tau);
	K_LP = f * exp(norminv( limit_delta/DF_f ) * sigma_deltaP * sqrt(tau) + 0.5 * sigma_deltaP^2 * tau);

	result = [K_atm, K_deltaC, K_deltaP, K_LC, K_LP];
end


function result = spot_pa(sigma_deltaC, sigma_deltaP, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv)

   global RR_delta_val;
	global limit_delta;
	
	e = 0.0001;
	E = 100 * S_0;
	
	 
	K_deltaP = fzero(@(K)spot_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, RR_delta_val), [e E]);
	Kmin = fmax(@(K)spot_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, RR_delta_val), [e E]);
	K_deltaC = fzero(@(K)spot_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, RR_delta_val), [Kmin K_deltaP]);
	K_atm = ATM_strike(S_0, f, sigma_atm, tau, ATM_conv, -1);
	K_LC = fzero(@(K)spot_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, limit_delta), [Kmin K_deltaP]); 
	K_LP = fzero(@(K)spot_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, limit_delta), [e E]);

	result = [K_atm, K_deltaC, K_deltaP, K_LC, K_LP];
end


function result = fwd_un(sigma_deltaC, sigma_deltaP, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv)

	global RR_delta_val;
	global limit_delta;

	K_deltaC = f * exp(-norminv(RR_delta_val) * sigma_deltaC * sqrt(tau) + 0.5 * sigma_deltaC^2 * tau);
	K_deltaP = f * exp(norminv(RR_delta_val) * sigma_deltaP * sqrt(tau) + 0.5 * sigma_deltaP^2 * tau);
	K_atm = ATM_strike(S_0, f, sigma_atm, tau, ATM_conv, 1);
	K_LC = f * exp(-norminv(limit_delta) * sigma_deltaC * sqrt(tau) + 0.5 * sigma_deltaC^2 * tau);
	K_LP = f * exp(norminv(limit_delta) * sigma_deltaP * sqrt(tau) + 0.5 * sigma_deltaP^2 * tau);

	result = [K_atm, K_deltaC, K_deltaP, K_LC, K_LP];
end


function result = fwd_pa(sigma_deltaC, sigma_deltaP, sigma_atm, S_0, f, tau, DF_f, DF_d, ATM_conv)

	global RR_delta_val;
	global limit_delta;
	
	e = 0.0001;
	E = 100 * S_0;

	K_deltaP = fzero(@(K)fwd_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, RR_delta_val), [e E]);
	Kmin = fmax(@(K)fwd_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, RR_delta_val), [e E]);
	K_deltaC = fzero(@(K)fwd_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, RR_delta_val), [Kmin K_deltaP]);
	K_atm = ATM_strike(S_0, f, sigma_atm, tau, ATM_conv, -1);
	K_LC = fzero(@(K)fwd_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, limit_delta), [Kmin K_deltaP]); 
	K_LP = fzero(@(K)fwd_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, limit_delta), [e E]);


	result = [K_atm, K_deltaC, K_deltaP, K_LC, K_LP];
end

% funkcja 12 ---------------------------------------------------------
% opis ATM_strike: funkcja wylicza strike K_atm dla roznych konwencji  
% input: S_0, f, sigma_atm, tau, 
% ATM_conv = {delta_neutral, atm_fwd, atm_spot, atm_value_neutral}, pa = 1 (unadjusted); -1 (premium adjusted)  
% output: K_atm

function K_atm = ATM_strike(S_0, f, sigma_atm, tau, ATM_conv, pa)
	if strcmp(ATM_conv, "delta_neutral") == 1
		K_atm =  f * exp((0.5 * pa * sigma_atm^2)*tau);
	elseif strcmp(ATM_conv, "atm_fwd")== 1
		K_atm = f;
	elseif strcmp(ATM_conv, "atm_spot") == 1 
		K_atm = S_0;
	elseif strcmp(ATM_conv, "atm_value_neutral") == 1
		K_atm = f;
	end
end

% funkcje 13-16 ------------------------------------------------------
% opis fwd_pa_call, fwd_pa_put, spot_pa_call, spot_pa_put: 
% funkcje pomocnicze wykorzystywane do obliczeñ spot_pa i fwd_pa
% input: K, f, sigma_deltaC, tau, DF_f, DF_d

function result = fwd_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, delta_val)
	result = K / f * (normcdf((log(f / K) - 0.5 * sigma_deltaC^2 * tau) / (sigma_deltaC * sqrt(tau)))) - delta_val;
endfunction

function result = fwd_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, delta_val)
	result =   K / f * (normcdf(-((log(f / K) - 0.5 * sigma_deltaP^2 * tau)) / (sigma_deltaP * sqrt(tau)))) - delta_val;
endfunction

function result = spot_pa_call(K, f, sigma_deltaC, tau, DF_f, DF_d, delta_val)
	result = DF_f * K / f * (normcdf((log(f / K) - 0.5 * sigma_deltaC ^ 2 * tau) / (sigma_deltaC * sqrt(tau)))) - delta_val;
endfunction

function result = spot_pa_put(K, f, sigma_deltaP, tau, DF_f, DF_d, delta_val)
	result =  DF_f * K / f * (normcdf(-((log(f / K) - 0.5 * sigma_deltaP^2 * tau)) / (sigma_deltaP * sqrt(tau)))) - delta_val;
endfunction

% funkcja 17 ------------------------------------------------
% opis ImpVol: zwraca wartosc zmiennosci implikowanej dla danego strika, czasu i wybranych metod interpolacyjnych (zmienna time_interp: typ interpolacji
% wzglêdem czasu jest zmienn¹ globalnym¹)
% input: 
%   strike - cena wykonania
%   start_date - data pocz¹tkowa
%   expire_date - data zapadalnosci
% output: wartosæ zmiennosci implikowanej

function result = ImpVol(strike, start_date, expire_date)
  
  vols = VolsInTime(strike);
  result = TimeInterpolate(start_date, expire_date, vols);

endfunction

% funkcja 18-----------------------------------------------------
% opis TimeInterpolate: dla danych na wejsciu zmiennosci dla tego samego strike'a w kolejnych momentach czasu oblicza wartosc zmiennosci implikowanej dla
% tego samego strike'a oraz czasu podanego na wejsciu
% input: 
%   start_date
%   expire_date
%   vols - wektor volatilities otrzymanych z krzywej volatility smile dla danego strike'a w poszczególnych momentach czasu w T_grid (siatce po czasie)
% output: wartosc zmiennosci implikowanej dla okresu miêdzy start_date a expire_date

function result = TimeInterpolate(start_date, expire_date, vols)

  global VOL_DCC;
  global volatility_grid;
  global T_grid;
  global smile_interp;
  global time_interp;
  
  tau = year_frac(start_date, expire_date, VOL_DCC);
  
  if strcmp(time_interp, "linear_on_variance") == 1
  % metoda interpolacji "linear_on_variance" polega na liniowej interpolacji funkcji t -> sigma(t)^2*t dla t pomiêdzy pierwszym a ostatnim wêz³em w T_grid
  % oraz przyjêciu sta³ej wartosci volatility dla czasu t poza siatka T_grid (równej wartosci w pierwszym i ostatnim wezle odpowiednio)
  
     index_T = 0;
	 
     for i = 1 : length(T_grid) 
       if (tau > T_grid(i))
         index_T = i;
       endif   
     endfor
	 
	 if (index_T == 0)
	   result = vols(1);
	 elseif (index_T == length(T_grid))
	   result = vols(length(T_grid));
	 else
	   TL = index_T;
	   TU = TL + 1;
	 
	   YL = vols(TL)^2*T_grid(TL);
	   YU = vols(TU)^2*T_grid(TU);
	 
	   if  (T_grid(TU) ~= T_grid(TL))       
         Y =  YL *  (T_grid(TU) -tau)/(T_grid(TU)-T_grid(TL)) +  YU * (tau - T_grid(TL))/(T_grid(TU)-T_grid(TL));  
       else
         Y = YL;
       endif

       result = sqrt(Y/tau);	   
	   
     endif
	
  else %domyslnie interpolacja liniowa - zwyk³a liniowa interpolacja funkcji t -> sigma(t)
  
     index_T = 0;
     for i = 1 : length(T_grid) 
       if (tau > T_grid(i))
         index_T = i;
       endif   
     endfor
	 if (index_T == 0)
	   TL = 1;
	   TU = 1;
	 elseif (index_T == length(T_grid))
	   TL = length(T_grid);
	   TU = length(T_grid);
	 else
	   TL = index_T;
	   TU = TL + 1;
	 endif
	 
	 YL = vols(TL);
	 YU = vols(TU);
	 
	 if  (T_grid(TU) ~= T_grid(TL))       
       Y =  YL *  (T_grid(TU) -tau)/(T_grid(TU)-T_grid(TL)) +  YU * (tau - T_grid(TL))/(T_grid(TU)-T_grid(TL));  
     else
       Y = YL;
     endif  
  
     result =  Y;	
	 
  endif
  
endfunction

% funkcja 19-------------------------------------------------
% opis VolsInTime: zwraca wektor volatilities otrzymanych z krzywej volatility smile dla danego strike'a w poszczególnych momentach czasu w T_grid 
%(siatce po czasie), przygotowuje input dla funkcji TimeInterpolate
%
% input: strike
% output: vols (wektor volatilities)

function vols = VolsInTime(strike)

  global volatility_grid;
  global K_grid;
  global T_grid;
  global m_points;
  global smile_interp;
  
  vols = zeros(1,length(T_grid));
  
  for t = 1:length(T_grid)
  
     index_K_L = 0;
     for i = 1 : columns(K_grid) 
       if (strike > K_grid(t,i))
        index_K_L = i;
       endif   
     endfor
     
     if (index_K_L == 0) 
         K_L = 1;
         K_R = 1;
     elseif (index_K_L == m_points)
	     K_L = m_points;
         K_R = m_points; 
     else
	     K_L = index_K_L;
         K_R = index_K_L+1;
     endif      
  
% linear interpolation between two strikes on the left and right of a given strike
% this is OK 
% strike interpolation is hiden in the construction of volatility_grid
  
     if (K_grid(t,K_R) ~= K_grid(t,K_L))
         Y = volatility_grid(t,K_L) * (K_grid(t,K_R)- strike)/(K_grid(t,K_R)-K_grid(t,K_L))...
         + volatility_grid(t,K_R) * (strike - K_grid(t,K_L))/(K_grid(t,K_R)-K_grid(t,K_L));
     else
         Y = volatility_grid(t,K_L);
     endif
     vols(t) = Y;
   endfor
  
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USMIECHY - RÓ¯NE RYSUNKI

% opis vol_surface: funkcja oblicza i rysuje powierzchniê zmiennosci implikowanej 
% input: no inputs ale obliczenia sa robione dla ustalonego smile_interp - interpolacja dla strików, 
% time_interp - interpolacja dla czasu jest liniowa
% output: powierzchnia zmiennosci implikowanej
 
 function   vol_surface( )
  global vol_fx_rate;
  global T_grid;
  global m_points;
  global t_points;
  global volatility_grid;
  global K_grid;
  global bandW_K;
  
  mm_points = m_points;
  T_plot = linspace(min(T_grid), max(T_grid), t_points);
  K_min = min(K_grid(:,1)); 
  K_max = max(K_grid(:,m_points));
   
  K_plot= linspace(K_min, K_max, mm_points);
  vol_surf = zeros(t_points, mm_points);

  for i = 1 : t_points 
  		tau = T_plot(i);
  		index_T = 0;
  		for j = 1 : length(T_grid) 
     		if (tau > T_grid(j))
        		index_T = j;
     		endif   
 		endfor
    
  		if (index_T >0 & index_T < length(T_grid))
			tL  = T_grid(index_T);
  			tU  = T_grid(index_T+1);
  			I_L = 0;
  			I_U = 0;    
  			for j = 1 : mm_points
       		if (K_plot(j)>=K_grid(index_T,1) & K_plot(j)<= K_grid(index_T,m_points))
       			if (I_L ==0)
       				I_LL = j;
       				I_L = 1;
       			endif	
       			ZL(j) =  volatility_grid(index_T, :)  * ker_W(K_plot(j), K_grid(index_T, :), bandW_K)'...
       			/sum(ker_W(K_plot(j), K_grid(index_T, :), bandW_K)); 
       			I_LU = j;
       		else
       			ZL(j) = NaN;	
       		endif
       		 
       		if (K_plot(j)>=K_grid(index_T+1,1) & K_plot(j)<= K_grid(index_T+1,m_points))
       			if (I_U ==0)
       				I_UL = j;
       				I_U = 1;
       			endif		
       			ZU(j) =  volatility_grid(index_T+1, :)  * ker_W(K_plot(j), K_grid(index_T+1, :), bandW_K)'...
       			/sum(ker_W(K_plot(j), K_grid(index_T+1, :), bandW_K));
       			I_UU = j; 
       		else
       			ZU(j) = NaN;
       		endif		  
  			endfor
  			vol_surf(i,:) = ZL *  (tU -tau)/(tU-tL) +  ZU * (tau - tL)/(tU-tL);
  			
  			I_L_min = round((I_UL - I_LL) * (tau - tL)/(tU-tL)) + I_LL;
  			for j = I_L_min : I_LL
  			   I_V_L = max(round((j - I_LL) * (tU-tL)/(tau - tL)) + I_LL,1);
  				vol_surf(i,j) = ZL(I_LL) *  (tU -tau)/(tU-tL) +  ZU(I_V_L) * (tau - tL)/(tU-tL);
  			endfor
  			
  			I_U_max = round((I_UU - I_LU) * (tau - tL)/(tU-tL)) + I_LU;
  			for j =  I_LU : I_U_max
  				I_V_U = min(round((j - I_LU) * (tU-tL)/(tau - tL)) + I_LU, mm_points);
  				vol_surf(i,j) = ZL(I_LU) *  (tU -tau)/(tU-tL) +  ZU(I_V_U) * (tau - tL)/(tU-tL);
  			endfor
  			  		
   	endif
    
  		if (index_T == 0)
  			K_L = K_grid(1,1); 
  			K_U = K_grid(1,m_points);
    		for j = 1 : mm_points
       		if (K_plot(j)>=K_L & K_plot(j)<= K_U)
       			Z(j) =  volatility_grid(1, :)  * ker_W(K_plot(j), K_grid(1, :), bandW_K)'...
       			/sum(ker_W(K_plot(j), K_grid(1, :), bandW_K)); 
       		else
       			Z(j) = NaN;
       		endif		  
  			endfor
    		vol_surf(i,:) = Z;
  		endif
  
  		if (index_T == length(T_grid))
  			K_L = K_grid(length(T_grid),1); 
  			K_U = K_grid(length(T_grid),m_points);
    		for j = 1 : mm_points
       		if (K_plot(j)>=K_L & K_plot(j)<= K_U)
       			Z(j) =  volatility_grid(length(T_grid), :)  * ker_W(K_plot(j), K_grid(length(T_grid), :), bandW_K)'...
       			/sum(ker_W(K_plot(j), K_grid(length(T_grid), :), bandW_K)); 
       		else
       			Z(j) = NaN;
       		endif		  
  			endfor
    		vol_surf(i,:) = Z;
  		endif
 endfor 
   
  
  mesh( K_plot/vol_fx_rate, T_plot, vol_surf);
  ylabel("maturity");
  xlabel("moneyness");
  zlabel("volatility");
 % print("volatility_surface.eps", "-deps");
 % close; 
  
endfunction
   

% opis smile_i: funkcja rysuje i zwraca krzywa zmiennosci implikowanej dla i-tego wiersza danych z fx_vol_data

function  smile_i(index)
  global vol_fx_rate;
  global volatility_grid;
  global K_grid;
  global fx_vol_data;
   
  
  strikes = fx_vol_data(index, 1:3);
  sigmas = fx_vol_data(index, 4:6);
  DFs = fx_vol_data(index, 7:8);
  
 
  strike_grid = K_grid(index,:);
  Y = volatility_grid(index,:);
  
  DF_f = DFs(1);
  DF_d = DFs(2);
  f = vol_fx_rate * (DF_f/DF_d);
  X = strike_grid / f;
   
  plot(X, Y,"b-", strikes/f, sigmas, "ro", "markersize", 3);
  xlabel("moneyness");
  ylabel("volatility");
#  plot_name = strcat("volatility_", FX_VOL{index,1},".eps");
#  print(plot_name, "-deps");
#  close;
endfunction


% opis smile: rysuje i zwraca krzywa zmiennosci implikowanej dla danej daty zapadalnosci  

function  smile(fixdate)
  global start_date;
  global vol_start_date;
  global vol_fx_rate;
  global DF_QUOT;
  global DF_BASE; 
  global VOL_DCC;
  global T_grid;
  global m_points;
  global volatility_grid;
  global K_grid;
  global smile_interp;
  global bandW_K;
  
    
  DF_f = DF(start_date, {fixdate}, DF_BASE);
  DF_d = DF(start_date, {fixdate}, DF_QUOT);
  f = vol_fx_rate * DF_d/DF_f;
  tau = year_frac(vol_start_date, {fixdate}, VOL_DCC);
  
  index_T = 0;
  for i = 1 : length(T_grid) 
     if (tau > T_grid(i))
        index_T = i;
     endif   
  endfor
    
  if (index_T >0 & index_T < length(T_grid))
     K_min = max(K_grid(index_T,1),K_grid(index_T+1,1)); 
     K_max = min(K_grid(index_T,m_points),K_grid(index_T+1,m_points));
     K_plot= linspace(K_min, K_max, m_points);
     tL  = T_grid(index_T);
     tU  =  T_grid(index_T+1);
       
     for i = 1 : m_points
       ZL(i) =  volatility_grid(index_T, :)  * ker_W(K_plot(i), K_grid(index_T, :), bandW_K)'/sum(ker_W(K_plot(i), K_grid(index_T, :), bandW_K)); 
       ZU(i) =  volatility_grid(index_T+1, :)  * ker_W(K_plot(i), K_grid(index_T+1, :), bandW_K)'/sum(ker_W(K_plot(i), K_grid(index_T+1, :), bandW_K));   
     endfor
    
    Z =  ZL *  (tU -tau)/(tU-tL) +  ZU * (tau - tL)/(tU-tL);
       
  endif  
  if (index_T == 0)
    K_plot = K_grid(1, :);
    Z = volatility_grid(1, :);
  endif
  
  if (index_T == length(T_grid))
    K_plot = K_grid(length(T_grid), :);
    Z = volatility_grid(length(T_grid), :);
  endif
    
  X = K_plot /f;
  
 % plot(X, Y);
 % figure;
  plot(X,Z);
  xlabel("moneyness");
  ylabel("volatility");
 % plot_name = strcat("volatility_", fixdate,".eps");
 % print(plot_name, "-deps");
 % close; 
  
endfunction



###################################################################################


% METODY INTERPOLACYJNE

% 1.METODA VANNA - VOLGA

% opis vanna_volga_interpolate: metoda interpolacji dla ceny wykonania (strike) 
% przy danych trzech punktach [sigma_deltaP, K_deltaP], [sigma_ATM, K_ATM], [sigma_deltaC, K_deltaC]
% input: strikes, sigmas, DFs, tau, strike
% output: interpolowana wartosc zmiennosci dla ustalonego strika

function result = vanna_volga_interpolate(strikes, sigmas, DFs, tau, strike)
  global vol_fx_rate;

  K1 = strikes(1);
  K2 = strikes(2);
  K3 = strikes(3);
  sigma_25P = sigmas(1);
  sigma_atm = sigmas(2);
  sigma_25C = sigmas(3);
  DF_f = DFs(1);
  DF_d = DFs(2);
  f = vol_fx_rate * (DF_f/DF_d);
  
  vanna_volga_result = vanna_volga(strikes, sigmas, DFs, tau, strike);
  
  result = ffzero(@(sigma)(vanna_volga_result - C(f, strike, sigma,  tau,   DF_d)));
endfunction
 

% opis vanna_volga: funkcja pomocnicza dla metoda vanna-volga interpolacji zmiennosci 
% dla danych strików (opis teoretyczny w raporcie), korzysta z funkcji pomocniczych vega i C
% input: strikes, sigmas, DFs, tau, K 
% output: wartosc opcji kupna

function result = vanna_volga(strikes, sigmas, DFs, tau, K)

	global vol_fx_rate;

	K1 = strikes(1);
	K2 = strikes(2);
	K3 = strikes(3);
	sigma_25P = sigmas(1);
	sigma_atm = sigmas(2);
	sigma_25C = sigmas(3);
	DF_f = DFs(1);
	DF_d = DFs(2);

	S_0 = vol_fx_rate;
	f = S_0 * (DF_f/DF_d);

	x1 = vega(S_0, f, tau, sigma_atm, K, DF_f, DF_d) / vega(S_0, f, tau, sigma_atm, K1, DF_f, DF_d) * log(K2/K) * log(K3/K) / (log(K2/K1) * log(K3/K1));
	
	x2 = vega(S_0, f, tau, sigma_atm, K, DF_f, DF_d) / vega(S_0, f, tau, sigma_atm, K2, DF_f, DF_d) * log(K1/K) * log(K3/K) / (log(K1/K2) * log(K3/K2));
	
	x3 = vega(S_0, f, tau, sigma_atm, K, DF_f, DF_d) / vega(S_0, f, tau, sigma_atm, K3, DF_f, DF_d) * log(K1/K) * log(K2/K) / (log(K1/K3) * log(K2/K3));
	
	a = C(f, K, sigma_atm, tau, DF_d);
	b = C(f, K1, sigma_25P, tau, DF_d);
	c = C(f, K1, sigma_atm, tau, DF_d);
	d = C(f, K3, sigma_25C, tau, DF_d);
	e = C(f, K3, sigma_atm, tau, DF_d);
	
	result = a + x1 * (b - c) + x3 * (d - e);
end



% opis vega: funkcja pomocnicza dla funkcji vanna_volga, oblicza wartosci vegi, vanny i volgi
% input: S_0, f, tau, sigma, K, DF_f, DF_d
% output: wartosæ vega

function result = vega(S_0, f, tau, sigma, K, DF_f, DF_d)

d1 = (log(f / K) + tau * 0.5 * sigma^2) / (sigma * sqrt(tau));
d2 = d1 - sigma * sqrt(tau);

vega = S_0 * DF_f * sqrt(tau) * exp(-0.5 * d1^2) / sqrt(2 * pi);
vanna = vega * d1 * d2 / sigma;
volga =  - vega * d2 / (S_0 * sigma * sqrt(tau));

result = vega;
end



% opis C: funkcja pomocnicza dla funkcji vanna_volga, 
% oblicza wartosæ opcji kupna dla danych parametrów poczatkowych: S_0, K, sigma, f, tau, DF_f, DF_d

function result = C(f, K, sigma, tau, DF_d)
	d1 = (log(f / K) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau));
	d2 = d1 - sigma * sqrt(tau);

	result =   DF_d * ( f * normcdf(d1) - K * normcdf(d2));
end


% WYBOR METODY INTERPOLACYJNEJ
% opis strike_interpolate: wybiera metode interpolacjna dla zmiennosci i cen wykonania


function result = strike_interpolate(strikes, sigmas, DFs, tau, strike, smile_interp)
  if strcmp(smile_interp, "vanna-volga") == 1
    result = vanna_volga_interpolate(strikes, sigmas, DFs, tau, strike);
  else
    result = kernel_interpolate(strikes, sigmas, strike);
  endif
endfunction



% 2.METODA J¥DROWA 

% opis kernel_interpolate: metoda interpolacji (domyslnie wykorzystywana dla p³aszczyzny czasu i zmiennosci), 
% korzysta z funkcji gaussian_kernel (opis teoretyczny w raporcie)

function ker_weight = ker_W(x, X,  band_w, alpha, kernel)
    if (nargin <4) 
       alpha = 1;
    end   
    u = abs((X - x) / band_w/alpha);
    u_ind = u>1;
    u(u_ind) = 1;
    ker_weight = 15*(1 - u.^2).^2/16;
   % ker_weight = exp(- u .^ 2 / 2)  ;    
endfunction

function result = kernel_interpolate(X, Y, x, band_w, alpha, kernel)
     if (nargin <4)
       band_w = 0.03; 
       alpha = 1;   
    end   
    if (nargin <6) 
       kernel=@gaussian_kernel;
    end   
    ks = kernel((X - x) /band_w/alpha);
    result = sum(Y .* ks) / sum(ks);
endfunction

function result = gaussian_kernel(u, a = 5)
    result = exp(- u .^ 2 / (2 * a ^ 2));
endfunction






% FUNKCJE POMOCNICZE DO OBLICZANIA MIEJSC ZEROWYCH

% opis fmax: oblicza argument funkcji f, dla którego wartosæ pochodnej wynosi zero
% input: funkcja f, punkt startowy x0
% output: argument funkcji f, dla którego wartosæ pochodnej wynosi zero

function result = fmax(f, x0)
    result = fzero(@(x)derivative(f, x), x0);
end



% opis derivative: oblicza pochodna numeryczna dla funkcji f i ustalonego punktu x
% input: funkcja f, punkt x
% output: wartosæ pochodnej numerycznej

function result = derivative(f, x)
    h = 0.00001;
    result = (f(x + h) - f(x - h)) / (2 * h);
end



% opis ffzero: oblicza miejsce zerowe dla danej funkcji f
% input: funkcja f
% output: miejsce zerowe

function result = ffzero(f)
	x0 = NaN;
	h = 0.1;
	while h > 0.0001 && isnan(x0)
		for x = h : h : 3.0
			if f(x) >= 0
				x0 = x;
				break
			end
		end
	h /= 2;
	end

	if f(x0) == 0
		result = x0;
		else
		result = fzero(f, [x0 100]);
	end
end

%-----------------------------------------
% funkcja zapisujaca dane z rynku VOLATILITY do ponownej inicjalizacji

function init_vol_data_gen()

global VOL_DCC;
global m_points;
global t_points;
global limit_delta;
global smile_interp;
global time_interp; 
global DF_QUOT;
global DF_BASE; 
global vol_fx_rate;
global vol_start_date;
global bandW_K; 
global bandW_T;
global alpha;
   
  fid = fopen("init_vol_data.m", "w");
  fprintf(fid, "VOL_DCC = \"%s\";\n",  VOL_DCC);
  fprintf(fid, "m_points = %d;\n",  m_points);
  fprintf(fid, "t_points = %d;\n",  t_points);
  fprintf(fid, "limit_delta = %f;\n",  limit_delta);
  n = length(DF_QUOT);  
  fprintf(fid, "DF_QUOT ={");
  for i=1:n-1
    fprintf(fid, "\"%s\", %f, \"%s\";\n",DF_QUOT{i,1},DF_QUOT{i,2},DF_QUOT{i,3});
  end
  fprintf(fid, "\"%s\", %f, \"%s\"};\n",DF_QUOT{n,1},DF_QUOT{n,2},DF_QUOT{n,3});  
  n = length(DF_BASE);  
  fprintf(fid, "DF_BASE ={");
  for i=1:n-1
    fprintf(fid, "\"%s\", %f, \"%s\";\n",DF_BASE{i,1},DF_BASE{i,2},DF_BASE{i,3});
  end
  fprintf(fid, "\"%s\", %f, \"%s\"};\n",DF_BASE{n,1},DF_BASE{n,2},DF_BASE{n,3});  
  
  fprintf(fid, "vol_start_date = \"%s\";\n", vol_start_date);
  fprintf(fid, "vol_fx_rate = %f;\n",  vol_fx_rate);
  fprintf(fid, "bandW_K = %f;\n",  bandW_K);
  fprintf(fid, "bandW_T = %f;\n",  bandW_T);
  fprintf(fid, "alpha = %f;\n",  alpha);
  fprintf(fid, "smile_interp = \"%s\";\n", smile_interp);
  fprintf(fid, "time_interp = \"%s\";\n\n", time_interp);
  fclose(fid);
  
endfunction