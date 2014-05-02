%volatility data input and verification
 
 
global VOL_DCC;
global VOL_BDA;
global VOL_EMA;

global vol_start_date;
global RR_delta_val;
global FX_rate;
global vol_fx_rate_type; 
global vol_fx_rate;

global FX_VOL;
global premium_curr;
global fx_vol_data;
global limit_delta;
global m_points;
global t_points;
global volatility_grid;
global K_grid;
global T_grid;
global smile_interp;
global time_interp;
global limits_table;
global bandW_K;
global bandW_T;
global alpha; 

global DSQ_Ask;
global DSQ_Bid;
global DSQ_Ave;

global DSB_Ask;
global DSB_Bid;
global DSB_Ave;

 
global DF_QUOT;
global DF_BASE; 
global BAS_CURR;
global QUO_CURR;
global ST_BAS_QUO; 
global FC_DOM;
global FC_FOR;

volatility_data_in; %wejscie

%implied_volatility;

 
implied_volatility_ap_new;

    if (strcmp(vol_fx_rate_type,'Mean'))
    	vol_fx_rate = sum(FX_rate)/2;
    elseif (strcmp(vol_fx_rate_type,'Ask'))
    	vol_fx_rate = FX_rate(2);
    else
    	vol_fx_rate = FX_rate(1);
    endif;


%wygenerowanie tablicy market_vol. 
market_vol_gen();

 