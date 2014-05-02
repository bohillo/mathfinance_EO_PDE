vol_start_date = "24-aug-2009";
 



VOL_DCC = "ACT/360";
VOL_BDA = "mfbd";
VOL_EMA = "+1";

vol_fx_rate_type = "Mean";
RR_delta_val = 0.25;
m_points = 100;
t_points = 100;
limit_delta = 0.050000;

FX_VOL = {"ON",30,0.75,1.3;
"1W",30,0.75,1.3;
"1M",28,0.75,1.3;
"2M",25.5,0.75,1.3;
"3M",24,0.725,1.3;
"6M",22,0.725,1.3;
"1Y",18.5,0.725,1.3;
"2Y",16.5,0.7,1.3;
"3Y",14,0.7,1.3;
"5Y",12,0.65,1.3};

DF_QUOT = DSQ_Ave;
DF_BASE = DSB_Ave;

bandW_K = 0.010000;
bandW_T = 0.030000;
alpha = 5.000000;
smile_interp = "vanna-volga";
time_interp = "linear";