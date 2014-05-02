%foreign data input and verification
 
global CURR_FOR;
 
global FC_FOR;
 
global FOR_DCC;
global FOR_BDA;
global FOR_EMA; 
global FOR_day_to_spot;
global pip_val;
global PDR;
global FX_rate;
global SWAP_POINTS;
global BAS_CURR;
global QUO_CURR;
global CURR_PAIR;
 


 
global DSF_Ask;
global DSF_Bid;
global DSF_Ave;

global DSQ_Ask;
global DSQ_Bid;
global DSQ_Ave;

global DSB_Ask;
global DSB_Bid;
global DSB_Ave;

foreign_data_in; %wejscie

 


 

%wygenerowanie tablic market_fx...
 
market_fx=market_fx_gen();


 
