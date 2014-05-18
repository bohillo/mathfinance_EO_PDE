
addpath("input");
DOM_data_in;
DOM_dis_factor;
FOR_data_in;
FOR_dis_factor;
exotic_options;
format long;
global DCC;
global Mx;
global Mt;
global dsigma;

DCC = 'ACT/ACT';
Mx = 2000;
Mt = 500;
dsigma = 0.001;

S0 = 100;

barrier = 120;
strike = 100;
MF = 1;
issue_date = '25-Aug-2009';
expire_date = '06-Nov-2009';
window_start_date = '27-Aug-2009';
window_end_date = '04-Nov-2009';
 
PPO = 0;
OSO = 0;
type = 'bid';

F_bid = S0*exp(0.025952 * year_frac(issue_date, expire_date, "ACT/ACT")) ;
F_ask = F_bid;

monitoring_dates = repmat(cellstr(expire_date), [1, day_diff(issue_date, expire_date, 'ACT')]);
%% monitoring days
MF = 1;
i = 1;
date = issue_date;


do 
  monitoring_dates(i) = date;
  date = day_shift(date, FC_DOM, 1);
  i += 1;
until day_diff(date, expire_date, 'ACT') <= 0;



#DM_out(F_bid,F_ask,barrier,strike, monitoring_dates, \
 #      issue_date,expire_date,PPO,OSO,type, "up", "call")

#DM_in(F_bid,F_ask,barrier,strike, monitoring_dates, \
 #     issue_date,expire_date,PPO,OSO,type, "up", "call")
Window_out(F_bid,F_ask,barrier,strike, [], \
       issue_date, window_start_date, window_end_date, \
	   expire_date,PPO,OSO,type, "up", "call")
uoamcall(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)(1)
call(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)(1)
