%domestic data input and verification
clear all;
global CURR_DOM;
 
global FC_DOM;
 
global start_date;

global CURR_FRAFUT;


global Depo_DCC;
global Depo_BDA;
global Depo_day_to_spot;
global Depo_EMA;
global Depo_rates;


global FRA_DCC;
global FRA_BDA;
global FRA_day_to_spot;
global FRA_EMA;
global FRA_rates;

global IRS_DCC;
global IRS_BDA;
global IRS_day_to_spot;
global IRS_EMA;
global IRS_cf;
global IRS_rates;


global interp_method;


global DSD_Bid;
global DSD_Ask;
global DSD_Ave;

%wejscie albo FRA - 1 albo Futures - 2

%domestic_data_in1; 


domestic_data_in1;
%domestic_data_in3;

%domestic_data_in4;

funkcje_kalendarzowe;

discount_functions_nonewline;

%wygenerowanie tablicy market_ir...
market_ir=market_ir_gen();
 


