% replaced DSD by DSQ - quote currency discount factor
% DSF by DSB - base currency discount factor
% for portfolio of options which prices are differences of simple options
% it is checked if the result is a positive number
% if not option value and Greecs are put to zero
% 27.10 correction in option diamput

format long;

function [result_forward] = forward(x_bid,x_ask,issue_date,expire_date)
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 DF_d_bid = DF(issue_date,{expire_date},DSQ_Bid)(1);
 DF_d_ask = DF(issue_date,{expire_date},DSQ_Ask)(1);
 DF_f_bid = DF(issue_date,{expire_date},DSB_Bid)(1);
 DF_f_ask = DF(issue_date,{expire_date},DSB_Ask)(1);
 result_forward(1) = x_bid*DF_f_ask/DF_d_bid;
 result_forward(2) = x_ask*DF_f_bid/DF_d_ask;
endfunction

function [result_DFd] = differ_DFd(DF_d,tau)
 result_DFd(1) = 0;
 result_DFd(2) = 0;
 result_DFd(3) = 0;
 result_DFd(4) = 0;
 result_DFd(5) = -DF_d*(log(DF_d))/tau;
 result_DFd(6) = 0;
endfunction

function [result_F] = differ_F(F,DF_d,DF_f,tau,sigma)
 result_F(1) = DF_f/DF_d;
 result_F(2) = 1;
 result_F(3) = 0;
 result_F(4) = 0;
 result_F(5) = F*(log(DF_d) - log(DF_f))/tau;
 result_F(6) = 0;
endfunction

function [result_d] = differ_d(F,DF_d,DF_f,tau,sigma,strike,eta)
 result_d(1) = (log(F/strike) + 0.5*eta*(sigma^2)*tau)/(sigma*sqrt(tau));
 result_d(2) = DF_f/(DF_d*F*sigma*sqrt(tau));
 result_d(3) = 1/(F*sigma*sqrt(tau));
 result_d(4) = -(DF_f^2)/((DF_d^2)*(F^2)*sigma*sqrt(tau));
 result_d(5) = -1/((F^2)*sigma*sqrt(tau));
 result_d(6) = (0.5*sigma*result_d(1) - ((log(DF_f) - log(DF_d))/tau + 0.5*eta*(sigma^2))*sqrt(tau))/(sigma*tau);
 result_d(7) = (eta*sigma*sqrt(tau) - result_d(1))/sigma;
endfunction

function [result_h] = differ_h(F,DF_d,DF_f,tau,sigma,barrier,strike,omega)
 result_h(1) = (log(((DF_f^2)*(barrier^2))/((DF_d^2)*F*strike)) + 0.5*omega*(sigma^2)*tau)/(sigma*sqrt(tau));
 result_h(2) = -DF_f/(DF_d*F*sigma*sqrt(tau));
 result_h(3) = -1/(F*sigma*sqrt(tau));
 result_h(4) = (DF_f^2)/((DF_d^2)*(F^2)*sigma*sqrt(tau));
 result_h(5) = 1/((F^2)*sigma*sqrt(tau));
 result_h(6) = (0.5*sigma*result_h(1) - ((log(DF_f) - log(DF_d))/tau + 0.5*omega*(sigma^2))*sqrt(tau))/(sigma*tau);
 result_h(7) = (omega*sigma*sqrt(tau) - result_h(1))/sigma;
endfunction

function [result_Nd] = differ_Nd(F,DF_d,DF_f,tau,sigma,strike,phi,eta)
 result_d = differ_d(F,DF_d,DF_f,tau,sigma,strike,eta);
 d = result_d(1);
 result_Nd(1) = normcdf(phi*d,0,1);
 result_Nd(2) = phi*DF_f*normpdf(phi*d,0,1)/(DF_d*F*sigma*sqrt(tau));
 result_Nd(3) = phi*normpdf(phi*d,0,1)/(F*sigma*sqrt(tau));
 result_Nd(4) = -phi*(DF_f^2)*normpdf(phi*d,0,1)*(1 + (phi^2)*d/(sigma*sqrt(tau)))/((DF_d^2)*(F^2)*sigma*sqrt(tau));
 result_Nd(5) = -phi*normpdf(phi*d,0,1)*(1 + (phi^2)*d/(sigma*sqrt(tau)))/((F^2)*sigma*sqrt(tau));
 result_Nd(6) = phi*normpdf(phi*d,0,1)*(0.5*sigma*d - ((log(DF_f) - log(DF_d))/tau + 0.5*eta*(sigma^2))*sqrt(tau))/(sigma*tau);
 result_Nd(7) = phi*normpdf(phi*d,0,1)*(eta*sigma*sqrt(tau) - d)/sigma;
endfunction

function [result_Nh] = differ_Nh(F,DF_d,DF_f,tau,sigma,barrier,strike,eta,omega)
 result_h = differ_h(F,DF_d,DF_f,tau,sigma,barrier,strike,omega);
 h = result_h(1);
 result_Nh(1) = normcdf(eta*h,0,1);
 result_Nh(2) = -eta*DF_f*normpdf(eta*h,0,1)/(DF_d*F*sigma*sqrt(tau));
 result_Nh(3) = -eta*normpdf(eta*h,0,1)/(F*sigma*sqrt(tau));
 result_Nh(4) = eta*(DF_f^2)*normpdf(eta*h,0,1)*(1 - (eta^2)*h/(sigma*sqrt(tau)))/((DF_d^2)*(F^2)*sigma*sqrt(tau));
 result_Nh(5) = eta*normpdf(eta*h,0,1)*(1 - (eta^2)*h/(sigma*sqrt(tau)))/((F^2)*sigma*sqrt(tau));
 result_Nh(6) = eta*normpdf(eta*h,0,1)*(0.5*sigma*h - ((log(DF_f) - log(DF_d))/tau + 0.5*omega*(sigma^2))*sqrt(tau))/(sigma*tau);
 result_Nh(7) = eta*normpdf(eta*h,0,1)*(omega*sigma*sqrt(tau) - h)/sigma;
endfunction

function [result_l] = differ_l(F,DF_d,DF_f,tau,sigma,barrier,omega)
 w = omega + 2*(log(DF_f) - log(DF_d))/((sigma^2)*tau);
 result_l(1) = ((DF_f*barrier)/(DF_d*F))^w;
 result_l(2) = -w*(barrier^w)*((DF_d*F/DF_f)^(-w - 1));
 result_l(3) = -w*((DF_f*barrier/DF_d)^w)*(F^(-w - 1));
 result_l(4) = w*(w + 1)*(barrier^w)*((DF_d*F/DF_f)^(-w - 2));
 result_l(5) = w*(w + 1)*((DF_f*barrier/DF_d)^w)*(F^(-w - 2));
 result_l(6) = 0;
 result_l(7) = -4*(log(DF_f) - log(DF_d))*((DF_f*barrier/DF_d*F)^w)*log(DF_f*barrier/DF_d*F)/((sigma^3)*tau);
endfunction

function [result_B1] = differ_B1(F,DF_d,DF_f,tau,sigma,strike,phi,eta)
 result_Nd = differ_Nd(F,DF_d,DF_f,tau,sigma,strike,phi,eta);
 result_B1 = phi*result_Nd;
endfunction

function [result_B2] = differ_B2(F,DF_d,DF_f,tau,sigma,barrier,strike,phi,eta,omega)
 result_Nh = differ_Nh(F,DF_d,DF_f,tau,sigma,barrier,strike,eta,omega);
 Nh = result_Nh(1);
 dNh_dx = result_Nh(2);
 dNh_dF = result_Nh(3);
 d2Nh_dx2 = result_Nh(4);
 d2Nh_dF2 = result_Nh(5);
 dNh_dt = result_Nh(6);
 dNh_dsigma = result_Nh(7);
 result_l = differ_l(F,DF_d,DF_f,tau,sigma,barrier,omega);
 l = result_l(1);
 dl_dx = result_l(2);
 dl_dF = result_l(3);
 d2l_dx2 = result_l(4);
 d2l_dF2 = result_l(5);
 dl_dt = result_l(6);
 dl_dsigma = result_l(7);
 result_B2(1) = phi*Nh*l;
 result_B2(2) = phi*(dNh_dx*l + Nh*dl_dx);
 result_B2(3) = phi*(dNh_dF*l + Nh*dl_dF);
 result_B2(4) = phi*(d2Nh_dx2*l + 2*dNh_dx*dl_dx + Nh*d2l_dx2);
 result_B2(5) = phi*(d2Nh_dF2*l + 2*dNh_dF*dl_dF + Nh*d2l_dF2);
 result_B2(6) = phi*(dNh_dt*l + Nh*dl_dt);
 result_B2(7) = phi*(dNh_dsigma*l + Nh*dl_dsigma);
endfunction

function [result_A1] = differ_A1(F,DF_d,DF_f,tau,sigma,strike,phi)
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1p = differ_B1(F,DF_d,DF_f,tau,sigma,strike,phi,1);
 B1p = result_B1p(1);
 dB1_dxp = result_B1p(2);
 dB1_dFp = result_B1p(3);
 d2B1_dx2p = result_B1p(4);
 d2B1_dF2p = result_B1p(5);
 dB1_dtp = result_B1p(6);
 dB1_dsigmap = result_B1p(7);
 result_B1d = differ_B1(F,DF_d,DF_f,tau,sigma,strike,phi,-1);
 B1d = result_B1d(1);
 dB1_dxd = result_B1d(2);
 dB1_dFd = result_B1d(3);
 d2B1_dx2d = result_B1d(4);
 d2B1_dF2d = result_B1d(5);
 dB1_dtd = result_B1d(6);
 dB1_dsigmad = result_B1d(7);
 result_A1(1) = F*B1p - strike*B1d;
 result_A1(2) = dF_dx*B1p + F*dB1_dxp - strike*dB1_dxd;
 result_A1(3) = dF_dF*B1p + F*dB1_dFp - strike*dB1_dFd;
 result_A1(4) = d2F_dx2*B1p + 2*dF_dx*dB1_dxp + F*d2B1_dx2p - strike*d2B1_dx2d;
 result_A1(5) = d2F_dF2*B1p + 2*dF_dF*dB1_dFp + F*d2B1_dF2p - strike*d2B1_dF2d;
 result_A1(6) = dF_dt*B1p + F*dB1_dtp - strike*dB1_dtd;
 result_A1(7) = dF_dsigma*B1p + F*dB1_dsigmap - strike*dB1_dsigmad;
endfunction

function [result_A2] = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,phi)
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1p = differ_B1(F,DF_d,DF_f,tau,sigma,barrier,phi,1);
 B1p = result_B1p(1);
 dB1_dxp = result_B1p(2);
 dB1_dFp = result_B1p(3);
 d2B1_dx2p = result_B1p(4);
 d2B1_dF2p = result_B1p(5);
 dB1_dtp = result_B1p(6);
 dB1_dsigmap = result_B1p(7);
 result_B1d = differ_B1(F,DF_d,DF_f,tau,sigma,barrier,phi,-1);
 B1d = result_B1d(1);
 dB1_dxd = result_B1d(2);
 dB1_dFd = result_B1d(3);
 d2B1_dx2d = result_B1d(4);
 d2B1_dF2d = result_B1d(5);
 dB1_dtd = result_B1d(6);
 dB1_dsigmad = result_B1d(7);
 result_A2(1) = F*B1p - strike*B1d;
 result_A2(2) = dF_dx*B1p + F*dB1_dxp - strike*dB1_dxd;
 result_A2(3) = dF_dF*B1p + F*dB1_dFp - strike*dB1_dFd;
 result_A2(4) = d2F_dx2*B1p + 2*dF_dx*dB1_dxp + F*d2B1_dx2p - strike*d2B1_dx2d;
 result_A2(5) = d2F_dF2*B1p + 2*dF_dF*dB1_dFp + F*d2B1_dF2p - strike*d2B1_dF2d;
 result_A2(6) = dF_dt*B1p + F*dB1_dtp - strike*dB1_dtd;
 result_A2(7) = dF_dsigma*B1p + F*dB1_dsigmap - strike*dB1_dsigmad;
endfunction

function [result_A3] = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,phi,eta)
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B2p = differ_B2(F,DF_d,DF_f,tau,sigma,strike,barrier,phi,eta,1);
 B2p = result_B2p(1);
 dB2_dxp = result_B2p(2);
 dB2_dFp = result_B2p(3);
 d2B2_dx2p = result_B2p(4);
 d2B2_dF2p = result_B2p(5);
 dB2_dtp = result_B2p(6);
 dB2_dsigmap = result_B2p(7);
 result_B2d = differ_B2(F,DF_d,DF_f,tau,sigma,strike,barrier,phi,eta,-1);
 B2d = result_B2d(1);
 dB2_dxd = result_B2d(2);
 dB2_dFd = result_B2d(3);
 d2B2_dx2d = result_B2d(4);
 d2B2_dF2d = result_B2d(5);
 dB2_dtd = result_B2d(6);
 dB2_dsigmad = result_B2d(7);
 result_A3(1) = F*B2p - strike*B2d;
 result_A3(2) = dF_dx*B2p + F*dB2_dxp - strike*dB2_dxd;
 result_A3(3) = dF_dF*B2p + F*dB2_dFp - strike*dB2_dFd;
 result_A3(4) = d2F_dx2*B2p + 2*dF_dx*dB2_dxp + F*d2B2_dx2p - strike*d2B2_dx2d;
 result_A3(5) = d2F_dF2*B2p + 2*dF_dF*dB2_dFp + F*d2B2_dF2p - strike*d2B2_dF2d;
 result_A3(6) = dF_dt*B2p + F*dB2_dtp - strike*dB2_dtd;
 result_A3(7) = dF_dsigma*B2p + F*dB2_dsigmap - strike*dB2_dsigmad;
endfunction

function [result_A4] = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,phi,eta)
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B2p = differ_B2(F,DF_d,DF_f,tau,sigma,barrier,barrier,phi,eta,1);
 B2p = result_B2p(1);
 dB2_dxp = result_B2p(2);
 dB2_dFp = result_B2p(3);
 d2B2_dx2p = result_B2p(4);
 d2B2_dF2p = result_B2p(5);
 dB2_dtp = result_B2p(6);
 dB2_dsigmap = result_B2p(7);
 result_B2d = differ_B2(F,DF_d,DF_f,tau,sigma,barrier,barrier,phi,eta,-1);
 B2d = result_B2d(1);
 dB2_dxd = result_B2d(2);
 dB2_dFd = result_B2d(3);
 d2B2_dx2d = result_B2d(4);
 d2B2_dF2d = result_B2d(5);
 dB2_dtd = result_B2d(6);
 dB2_dsigmad = result_B2d(7);
 result_A4(1) = F*B2p - strike*B2d;
 result_A4(2) = dF_dx*B2p + F*dB2_dxp - strike*dB2_dxd;
 result_A4(3) = dF_dF*B2p + F*dB2_dFp - strike*dB2_dFd;
 result_A4(4) = d2F_dx2*B2p + 2*dF_dx*dB2_dxp + F*d2B2_dx2p - strike*d2B2_dx2d;
 result_A4(5) = d2F_dF2*B2p + 2*dF_dF*dB2_dFp + F*d2B2_dF2p - strike*d2B2_dF2d;
 result_A4(6) = dF_dt*B2p + F*dB2_dtp - strike*dB2_dtd;
 result_A4(7) = dF_dsigma*B2p + F*dB2_dsigmap - strike*dB2_dsigmad;
endfunction

function [result_call] = call(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1p = differ_B1(F,DF_d,DF_f,tau,sigma,strike,1,1);
 B1p = result_B1p(1);
 dB1_dxp = result_B1p(2);
 dB1_dFp = result_B1p(3);
 d2B1_dx2p = result_B1p(4);
 d2B1_dF2p = result_B1p(5);
 dB1_dtp = result_B1p(6);
 dB1_dsigmap = result_B1p(7);
 result_B1d = differ_B1(F,DF_d,DF_f,tau,sigma,strike,1,-1);
 B1d = result_B1d(1);
 dB1_dxd = result_B1d(2);
 dB1_dFd = result_B1d(3);
 d2B1_dx2d = result_B1d(4);
 d2B1_dF2d = result_B1d(5);
 dB1_dtd = result_B1d(6);
 dB1_dsigmad = result_B1d(7);
 result_call(1) = DF_dG*(F*B1p - strike*B1d);
 result_call(2) = DF_dG*(dF_dx*B1p + F*dB1_dxp - strike*dB1_dxd);
 result_call(3) = DF_dG*(dF_dF*B1p + F*dB1_dFp - strike*dB1_dFd);
 result_call(4) = DF_dG*(d2F_dx2*B1p + 2*dF_dx*dB1_dxp + F*d2B1_dx2p - strike*d2B1_dx2d);
 result_call(5) = DF_dG*(d2F_dF2*B1p + 2*dF_dF*dB1_dFp + F*d2B1_dF2p - strike*d2B1_dF2d);
 result_call(6) = dDFdG_dt*(F*B1p - strike*B1d) + DF_dG*(dF_dt*B1p + F*dB1_dtp - strike*dB1_dtd);
 result_call(7) = DF_dG*(dF_dsigma*B1p + F*dB1_dsigmap - strike*dB1_dsigmad);
 if (result_call(1)<=0)
    result_call = zeros(1,7);
 endif 
endfunction

function [result_put] = put(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1p = differ_B1(F,DF_d,DF_f,tau,sigma,strike,-1,1);
 B1p = result_B1p(1);
 dB1_dxp = result_B1p(2);
 dB1_dFp = result_B1p(3);
 d2B1_dx2p = result_B1p(4);
 d2B1_dF2p = result_B1p(5);
 dB1_dtp = result_B1p(6);
 dB1_dsigmap = result_B1p(7);
 result_B1d = differ_B1(F,DF_d,DF_f,tau,sigma,strike,-1,-1);
 B1d = result_B1d(1);
 dB1_dxd = result_B1d(2);
 dB1_dFd = result_B1d(3);
 d2B1_dx2d = result_B1d(4);
 d2B1_dF2d = result_B1d(5);
 dB1_dtd = result_B1d(6);
 dB1_dsigmad = result_B1d(7);
 result_put(1) = DF_dG*(F*B1p - strike*B1d);
 result_put(2) = DF_dG*(dF_dx*B1p + F*dB1_dxp - strike*dB1_dxd);
 result_put(3) = DF_dG*(dF_dF*B1p + F*dB1_dFp - strike*dB1_dFd);
 result_put(4) = DF_dG*(d2F_dx2*B1p + 2*dF_dx*dB1_dxp + F*d2B1_dx2p - strike*d2B1_dx2d);
 result_put(5) = DF_dG*(d2F_dF2*B1p + 2*dF_dF*dB1_dFp + F*d2B1_dF2p - strike*d2B1_dF2d);
 result_put(6) = dDFdG_dt*(F*B1p - strike*B1d) + DF_dG*(dF_dt*B1p + F*dB1_dtp - strike*dB1_dtd);
 result_put(7) = DF_dG*(dF_dsigma*B1p + F*dB1_dsigmap - strike*dB1_dsigmad);
 if (result_put(1)<=0)
    result_put = zeros(1,7);
 endif 
endfunction

function [result_riskrev] = riskrev(F_bid,F_ask,strike_1,strike_2,issue_date,expire_date,PPO,OSO,type)
 if (type == "bid")
   opp_type = "ask";
  else
   opp_type = "bid";
 endif
 result_riskrev = call(F_bid,F_ask,strike_2,issue_date,expire_date,PPO,OSO,type) - put(F_bid,F_ask,strike_1,issue_date,expire_date, PPO,OSO,opp_type);
 if (result_riskrev(1) <= 0)
    result_riskrev = zeros(1,7);
 endif    
endfunction

function [result_straddle] = straddle(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 result_straddle = call(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type) + put(F_bid,F_ask,strike,issue_date,expire_date, PPO,OSO,type);
 if (result_straddle(1)<=0)
    result_straddle = zeros(1,7);
 endif 
endfunction

function [result_strangle] = strangle(F_bid,F_ask,strike_1,strike_2,issue_date,expire_date,PPO,OSO,type)
 result_strangle = call(F_bid,F_ask,strike_2,issue_date,expire_date,PPO,OSO,type) + put(F_bid,F_ask,strike_1,issue_date,expire_date,PPO,OSO,type);
 if (result_strangle(1)<=0)
    result_strangle = zeros(1,7);
 endif 
endfunction

function [result_butterfly] = butterfly(F_bid,F_ask,strike_1,strike_2,issue_date,expire_date,PPO,OSO,type)
 if (type == "bid")
   opp_type = "ask";
  else
   opp_type = "bid";
 endif
 result_butterfly = call(F_bid,F_ask,strike_1,issue_date,expire_date,PPO,OSO,type) + call(F_bid,F_ask,strike_2,issue_date,expire_date,PPO,OSO,type) - 2*call(F_bid,F_ask,0.5*(strike_1 + strike_2),issue_date,expire_date, PPO,OSO,opp_type);
 if (result_butterfly(1) <= 0)
    result_butterfly = zeros(1,7);
 endif   
endfunction

function [result_seagull] = seagull(F_bid,F_ask,strike_1,strike_2,strike_3,issue_date,expire_date,PPO,OSO,type)
 if (type == "bid")
   opp_type = "ask";
  else
   opp_type = "bid";
 endif
 result_seagull = call(F_bid,F_ask,strike_2,issue_date,expire_date,PPO,OSO,type) - call(F_bid,F_ask,strike_3,issue_date,   expire_date,PPO,OSO,opp_type) - put(F_bid,F_ask,strike_1,issue_date,expire_date,PPO,OSO,opp_type);
 if (result_seagull(1) <= 0)
    result_seagull = zeros(1,7);
 endif
endfunction

function [result_concall] = concall(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1 = differ_B1(F,DF_d,DF_f,tau,sigma,strike,1,-1);
 B1 = result_B1(1);
 dB1_dx = result_B1(2);
 dB1_dF = result_B1(3);
 d2B1_dx2 = result_B1(4);
 d2B1_dF2 = result_B1(5);
 dB1_dt = result_B1(6);
 dB1_dsigma = result_B1(7);
 result_concall(1) = DF_dG*B1;
 result_concall(2) = DF_dG*dB1_dx;
 result_concall(3) = DF_dG*dB1_dF;
 result_concall(4) = DF_dG*d2B1_dx2;
 result_concall(5) = DF_dG*d2B1_dF2;
 result_concall(6) = dDFdG_dt*B1 + DF_dG*dB1_dt;
 result_concall(7) = DF_dG*dB1_dsigma;
 if (result_concall(1)<=0)
    result_concall = zeros(1,7);
 endif 
endfunction

function [result_conput] = conput(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1 = differ_B1(F,DF_d,DF_f,tau,sigma,strike,-1,-1);
 B1 = result_B1(1);
 dB1_dx = result_B1(2);
 dB1_dF = result_B1(3);
 d2B1_dx2 = result_B1(4);
 d2B1_dF2 = result_B1(5);
 dB1_dt = result_B1(6);
 dB1_dsigma = result_B1(7);
 result_conput(1) = -DF_dG*B1;
 result_conput(2) = -DF_dG*dB1_dx;
 result_conput(3) = -DF_dG*dB1_dF;
 result_conput(4) = -DF_dG*d2B1_dx2;
 result_conput(5) = -DF_dG*d2B1_dF2;
 result_conput(6) = -dDFdG_dt*B1 - DF_dG*dB1_dt;
 result_conput(7) = -DF_dG*dB1_dsigma;
 if (result_conput(1)<=0)
    result_conput = zeros(1,7);
 endif 
endfunction

function [result_aoncall] = aoncall(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1 = differ_B1(F,DF_d,DF_f,tau,sigma,strike,1,1);
 B1 = result_B1(1);
 dB1_dx = result_B1(2);
 dB1_dF = result_B1(3);
 d2B1_dx2 = result_B1(4);
 d2B1_dF2 = result_B1(5);
 dB1_dt = result_B1(6);
 dB1_dsigma = result_B1(7);
 result_aoncall(1) = DF_dG*F*B1;
 result_aoncall(2) = DF_dG*(dF_dx*B1 + F*dB1_dx);
 result_aoncall(3) = DF_dG*(dF_dF*B1 + F*dB1_dF);
 result_aoncall(4) = DF_dG*(d2F_dx2*B1 + 2*dF_dx*dB1_dx + F*d2B1_dx2);
 result_aoncall(5) = DF_dG*(d2F_dF2*B1 + 2*dF_dF*dB1_dF + F*d2B1_dF2);
 result_aoncall(6) = dDFdG_dt*F*B1 + DF_dG*(dF_dt*B1 + F*dB1_dt);
 result_aoncall(7) = DF_dG*(dF_dsigma*B1 + F*dB1_dsigma);
 if (result_aoncall(1)<=0)
    result_aoncall = zeros(1,7);
 endif 
endfunction

function [result_aonput] = aonput(F_bid,F_ask,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 result_F = differ_F(F,DF_d,DF_f,tau,sigma);
 dF_dx = result_F(1);
 dF_dF = result_F(2);
 d2F_dx2 = result_F(3);
 d2F_dF2 = result_F(4);
 dF_dt = result_F(5);
 dF_dsigma = result_F(6);
 result_B1 = differ_B1(F,DF_d,DF_f,tau,sigma,strike,-1,1);
 B1 = result_B1(1);
 dB1_dx = result_B1(2);
 dB1_dF = result_B1(3);
 d2B1_dx2 = result_B1(4);
 d2B1_dF2 = result_B1(5);
 dB1_dt = result_B1(6);
 dB1_dsigma = result_B1(7);
 result_aonput(1) = -DF_dG*F*B1;
 result_aonput(2) = -DF_dG*(dF_dx*B1 + F*dB1_dx);
 result_aonput(3) = -DF_dG*(dF_dF*B1 + F*dB1_dF);
 result_aonput(4) = -DF_dG*(d2F_dx2*B1 + 2*dF_dx*dB1_dx + F*d2B1_dx2);
 result_aonput(5) = -DF_dG*(d2F_dF2*B1 + 2*dF_dF*dB1_dF + F*d2B1_dF2);
 result_aonput(6) = -dDFdG_dt*F*B1 - DF_dG*(dF_dt*B1 + F*dB1_dt);
 result_aonput(7) = -DF_dG*(dF_dsigma*B1 + F*dB1_dsigma);
 if (result_aonput(1)<=0)
    result_aonput = zeros(1,7);
 endif 
endfunction

function [result_doeucall] = doeucall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 if (barrier < strike)
   result_doeucall = call(F_bid,F_ask,strike,issue_date,expire_date, PPO,OSO,type);
  else
   result_doeucall = call(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,type) + (barrier - strike)*concall(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,type);
   
 endif
 if (result_doeucall(1)<=0)
    result_doeucall = zeros(1,7);
 endif 
  
endfunction

function [result_uieucall] = uieucall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 result_uieucall = doeucall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type);
endfunction

function [result_uoeucall] = uoeucall(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 if (type == "bid")
   opp_type = "ask";
  else
   opp_type = "bid";
 endif
 if (barrier >= strike)
   result_uoeucall = call(F_bid,F_ask,strike,issue_date,expire_date, PPO,OSO,type) - call(F_bid,F_ask,barrier,issue_date, expire_date, PPO,OSO,opp_type) + (strike - barrier)*concall(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,type);
  else
   result_uoeucall(1) = 0;
   result_uoeucall(2) = 0;
   result_uoeucall(3) = 0;
   result_uoeucall(4) = 0;
   result_uoeucall(5) = 0;
   result_uoeucall(6) = 0;
   result_uoeucall(7) = 0;
 endif
 if (result_uoeucall(1)<=0)
    result_uoeucall = zeros(1,7);
 endif 
endfunction

function [result_dieucall] = dieucall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 result_dieucall = uoeucall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type);
endfunction

function [result_doeuput] = doeuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 if (type == "bid")
   opp_type = "ask";
  else
   opp_type = "bid";
 endif
 if (barrier < strike)
   result_doeuput = put(F_bid,F_ask,strike,issue_date,expire_date, PPO,OSO,type) - put(F_bid,F_ask,barrier,issue_date,expire_date, OSO,opp_type) - (strike - barrier)*conput(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,opp_type);
  else
   result_doeuput(1) = 0;
   result_doeuput(2) = 0;
   result_doeuput(3) = 0;
   result_doeuput(4) = 0;
   result_doeuput(5) = 0;
   result_doeuput(6) = 0;
   result_doeuput(7) = 0;
 endif
 if (result_doeuput(1)<=0)
    result_doeuput = zeros(1,7);
 endif 
endfunction

function [result_uieuput] = uieuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 result_uieuput = doeuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type);
endfunction

function [result_uoeuput] = uoeuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 if (barrier < strike)
   result_uoeuput = put(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,type) + (strike - barrier)*conput(F_bid,F_ask,barrier,issue_date,expire_date, PPO,OSO,type);
  else
   result_uoeuput = put(F_bid,F_ask,strike,issue_date,expire_date, PPO,OSO,type);
 endif
 if (result_uoeuput(1)<=0)
    result_uoeuput = zeros(1,7);
 endif 
endfunction

function [result_dieuput] = dieuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 result_dieuput = uoeuput(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type);
endfunction

function [result_doamcall] = doamcall(F_bid,F_ask,barrier,strike,issue_date,expire_date, PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,1,1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_doamcall(1) = DF_dG*(A1 - A3);
   result_doamcall(2) = dDFdG_dx*(A1 - A3) + DF_dG*(dA1_dx - dA3_dx);
   result_doamcall(3) = dDFdG_dF*(A1 - A3) + DF_dG*(dA1_dF - dA3_dF);
   result_doamcall(4) = d2DFdG_dx2*(A1 - A3) + 2*dDFdG_dx*(dA1_dx - dA3_dx) + DF_dG*(d2A1_dx2 - d2A3_dx2);
   result_doamcall(5) = d2DFdG_dF2*(A1 - A3) + 2*dDFdG_dF*(dA1_dF - dA3_dF) + DF_dG*(d2A1_dF2 - d2A3_dF2);
   result_doamcall(6) = dDFdG_dt*(A1 - A3) + DF_dG*(dA1_dt - dA3_dt);
   result_doamcall(7) = dDFdG_dsigma*(A1 - A3) + DF_dG*(dA1_dsigma - dA3_dsigma);
  else
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,1,1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_doamcall(1) = DF_dG*(A2 - A4);
   result_doamcall(2) = dDFdG_dx*(A2 - A4) + DF_dG*(dA2_dx - dA4_dx);
   result_doamcall(3) = dDFdG_dF*(A2 - A4) + DF_dG*(dA2_dF - dA4_dF);
   result_doamcall(4) = d2DFdG_dx2*(A2 - A4) + 2*dDFdG_dx*(dA2_dx - dA4_dx) + DF_dG*(d2A2_dx2 - d2A4_dx2);
   result_doamcall(5) = d2DFdG_dF2*(A2 - A4) + 2*dDFdG_dF*(dA2_dF - dA4_dF) + DF_dG*(d2A2_dF2 - d2A4_dF2);
   result_doamcall(6) = dDFdG_dt*(A2 - A4) + DF_dG*(dA2_dt - dA4_dt);
   result_doamcall(7) = dDFdG_dsigma*(A2 - A4) + DF_dG*(dA2_dsigma - dA4_dsigma);
 endif
 if (result_doamcall(1)<=0)
    result_doamcall = zeros(1,7);
 endif    
endfunction

function [result_uiamcall] = uiamcall(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_uiamcall(1) = DF_dG*A1;
   result_uiamcall(2) = dDFdG_dx*A1 + DF_dG*dA1_dx;
   result_uiamcall(3) = dDFdG_dF*A1 + DF_dG*dA1_dF;
   result_uiamcall(4) = d2DFdG_dx2*A1 + 2*dDFdG_dx*dA1_dx + DF_dG*d2A1_dx2;
   result_uiamcall(5) = d2DFdG_dF2*A1 + 2*dDFdG_dF*dA1_dF + DF_dG*d2A1_dF2;
   result_uiamcall(6) = dDFdG_dt*A1 + DF_dG*dA1_dt;
   result_uiamcall(7) = dDFdG_dsigma*A1 + DF_dG*dA1_dsigma;
  else
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,1,-1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,1,-1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_uiamcall(1) = DF_dG*(A2 - A3 + A4);
   result_uiamcall(2) = dDFdG_dx*(A2 - A3 + A4) + DF_dG*(dA2_dx - dA3_dx + dA4_dx);
   result_uiamcall(3) = dDFdG_dF*(A2 - A3 + A4) + DF_dG*(dA2_dF - dA3_dF + dA4_dF);
   result_uiamcall(4) = d2DFdG_dx2*(A2 - A3 + A4) + 2*dDFdG_dx*(dA2_dx - dA3_dx + dA4_dx) + DF_dG*(d2A2_dx2 - d2A3_dx2 + d2A4_dx2);
   result_uiamcall(5) = d2DFdG_dF2*(A2 - A3 + A4) + 2*dDFdG_dF*(dA2_dF - dA3_dF + dA4_dF) + DF_dG*(d2A2_dF2 - d2A3_dF2 + d2A4_dF2);
   result_uiamcall(6) = dDFdG_dt*(A2 - A3 + A4) + DF_dG*(dA2_dt - dA3_dt + dA4_dt);
   result_uiamcall(7) = dDFdG_dsigma*(A2 - A3 + A4) + DF_dG*(dA3_dsigma - dA3_dsigma + dA4_dsigma);
 endif
 if (result_uiamcall(1)<=0)
    result_uiamcall = zeros(1,7);
 endif 
endfunction

function [result_uoamcall] = uoamcall(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_uoamcall(1) = 0;
   result_uoamcall(2) = 0;
   result_uoamcall(3) = 0;
   result_uoamcall(4) = 0;
   result_uoamcall(5) = 0;
   result_uoamcall(6) = 0;
   result_uoamcall(7) = 0;
  else
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,1,-1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,1,-1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_uoamcall(1) = DF_dG*(A1 - A2 + A3 - A4);
   result_uoamcall(2) = dDFdG_dx*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dx - dA2_dx + dA3_dx - dA4_dx);
   result_uoamcall(3) = dDFdG_dF*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dF - dA2_dF + dA3_dF - dA4_dF);
   result_uoamcall(4) = d2DFdG_dx2*(A1 - A2 + A3 - A4) + 2*dDFdG_dx*(dA1_dx - dA2_dx + dA3_dx - dA4_dx) + DF_dG*(d2A1_dx2 - d2A2_dx2 + d2A3_dx2 - d2A4_dx2);
   result_uoamcall(5) = d2DFdG_dF2*(A1 - A2 + A3 - A4) + 2*dDFdG_dF*(dA1_dF - dA2_dF + dA3_dF - dA4_dF) + DF_dG*(d2A1_dF2 - d2A2_dF2 + d2A3_dF2 - d2A4_dF2);
   result_uoamcall(6) = dDFdG_dt*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dt - dA2_dt + dA3_dt - dA4_dt);
   result_uoamcall(7) = dDFdG_dsigma*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dsigma - dA2_dsigma + dA3_dsigma - dA4_dsigma);
 endif
 if (result_uoamcall(1)<=0)
    result_uoamcall = zeros(1,7);
 endif 
endfunction

function [result_diamcall] = diamcall(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,1,1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_diamcall(1) = DF_dG*A3;
   result_diamcall(2) = dDFdG_dx*A3 + DF_dG*dA3_dx;
   result_diamcall(3) = dDFdG_dF*A3 + DF_dG*dA3_dF;
   result_diamcall(4) = d2DFdG_dx2*A3 + 2*dDFdG_dx*dA3_dx + DF_dG*d2A3_dx2;
   result_diamcall(5) = d2DFdG_dF2*A3 + 2*dDFdG_dF*dA3_dF + DF_dG*d2A3_dF2;
   result_diamcall(6) = dDFdG_dt*A3 + DF_dG*dA3_dt;
   result_diamcall(7) = dDFdG_dsigma*A3 + DF_dG*dA3_dsigma;
  else
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,1,1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_diamcall(1) = DF_dG*(A1 - A2 + A4);
   result_diamcall(2) = dDFdG_dx*(A1 - A2 + A4) + DF_dG*(dA1_dx - dA2_dx + dA4_dx);
   result_diamcall(3) = dDFdG_dF*(A1 - A2 + A4) + DF_dG*(dA1_dF - dA2_dF + dA4_dF);
   result_diamcall(4) = d2DFdG_dx2*(A1 - A2 + A4) + 2*dDFdG_dx*(dA1_dx - dA2_dx + dA4_dx) + DF_dG*(d2A1_dx2 - d2A2_dx2 + d2A4_dx2);
   result_diamcall(5) = d2DFdG_dF2*(A1 - A2 + A4) + 2*dDFdG_dF*(dA1_dF - dA2_dF + dA4_dF) + DF_dG*(d2A1_dF2 - d2A2_dF2 + d2A4_dF2);
   result_diamcall(6) = dDFdG_dt*(A1 - A2 + A4) + DF_dG*(dA1_dt - dA2_dt + dA4_dt);
   result_diamcall(7) = dDFdG_dsigma*(A1 - A2 + A4) + DF_dG*(dA1_dsigma - dA2_dsigma + dA4_dsigma);
 endif
 if (result_diamcall(1)<=0)
    result_diamcall = zeros(1,7);
 endif 
endfunction

function [result_doamput] = doamput(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,-1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,-1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_doamput(1) = DF_dG*(A1 - A2 + A3 - A4);
   result_doamput(2) = dDFdG_dx*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dx - dA2_dx + dA3_dx - dA4_dx);
   result_doamput(3) = dDFdG_dF*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dF - dA2_dF + dA3_dF - dA4_dF);
   result_doamput(4) = d2DFdG_dx2*(A1 - A2 + A3 - A4) + 2*dDFdG_dx*(dA1_dx - dA2_dx + dA3_dx - dA4_dx) + DF_dG*(d2A1_dx2 - d2A2_dx2 + d2A3_dx2 - d2A4_dx2);
   result_doamput(5) = d2DFdG_dF2*(A1 - A2 + A3 - A4) + 2*dDFdG_dF*(dA1_dF - dA2_dF + dA3_dF - dA4_dF) + DF_dG*(d2A1_dF2 - d2A2_dF2 + d2A3_dF2 - d2A4_dF2);
   result_doamput(6) = dDFdG_dt*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dt - dA2_dt + dA3_dt - dA4_dt);
   result_doamput(7) = dDFdG_dsigma*(A1 - A2 + A3 - A4) + DF_dG*(dA1_dsigma - dA2_dsigma + dA3_dsigma - dA4_dsigma); 
  else
   result_doamput(1) = 0;
   result_doamput(2) = 0;
   result_doamput(3) = 0;
   result_doamput(4) = 0;
   result_doamput(5) = 0;
   result_doamput(6) = 0;
   result_doamput(7) = 0;
 endif
 if (result_doamput(1)<=0)
    result_doamput = zeros(1,7);
 endif 
endfunction

function [result_uiamput] = uiamput(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,-1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,-1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,-1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_uiamput(1) = DF_dG*(A1 - A2 + A4);
   result_uiamput(2) = dDFdG_dx*(A1 - A2 + A4) + DF_dG*(dA1_dx - dA2_dx + dA4_dx);
   result_uiamput(3) = dDFdG_dF*(A1 - A2 + A4) + DF_dG*(dA1_dF - dA2_dF + dA4_dF);
   result_uiamput(4) = d2DFdG_dx2*(A1 - A2 + A4) + 2*dDFdG_dx*(dA1_dx - dA2_dx + dA4_dx) + DF_dG*(d2A1_dx2 - d2A2_dx2 + d2A4_dx2);
   result_uiamput(5) = d2DFdG_dF2*(A1 - A2 + A4) + 2*dDFdG_dF*(dA1_dF - dA2_dF + dA4_dF) + DF_dG*(d2A1_dF2 - d2A2_dF2 + d2A4_dF2);
   result_uiamput(6) = dDFdG_dt*(A1 - A2 + A4) + DF_dG*(dA1_dt - dA2_dt + dA4_dt);
   result_uiamput(7) = dDFdG_dsigma*(A1 - A2 + A4) + DF_dG*(dA1_dsigma - dA2_dsigma + dA4_dsigma);
  else
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,-1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_uiamput(1) = DF_dG*A3;
   result_uiamput(2) = dDFdG_dx*A3 + DF_dG*dA3_dx;
   result_uiamput(3) = dDFdG_dF*A3 + DF_dG*dA3_dF;
   result_uiamput(4) = d2DFdG_dx2*A3 + 2*dDFdG_dx*dA3_dx + DF_dG*d2A3_dx2;
   result_uiamput(5) = d2DFdG_dF2*A3 + 2*dDFdG_dF*dA3_dF + DF_dG*d2A3_dF2;
   result_uiamput(6) = dDFdG_dt*A3 + DF_dG*dA3_dt;
   result_uiamput(7) = dDFdG_dsigma*A3 + DF_dG*dA3_dsigma;
 endif
 if (result_uiamput(1)<=0)
    result_uiamput = zeros(1,7);
 endif 
endfunction

function [result_uoamput] = uoamput(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,-1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,-1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_uoamput(1) = DF_dG*(A2 - A4);
   result_uoamput(2) = dDFdG_dx*(A2 - A4) + DF_dG*(dA2_dx - dA4_dx);
   result_uoamput(3) = dDFdG_dF*(A2 - A4) + DF_dG*(dA2_dF - dA4_dF);
   result_uoamput(4) = d2DFdG_dx2*(A2 - A4) + 2*dDFdG_dx*(dA2_dx - dA4_dx) + DF_dG*(d2A2_dx2 - d2A4_dx2);
   result_uoamput(5) = d2DFdG_dF2*(A2 - A4) + 2*dDFdG_dF*(dA2_dF - dA4_dF) + DF_dG*(d2A2_dF2 - d2A4_dF2);
   result_uoamput(6) = dDFdG_dt*(A2 - A4) + DF_dG*(dA2_dt - dA4_dt);
   result_uoamput(7) = dDFdG_dsigma*(A2 - A4) + DF_dG*(dA2_dsigma - dA4_dsigma);
  else
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,-1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,-1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_uoamput(1) = DF_dG*(A1 - A3);
   result_uoamput(2) = dDFdG_dx*(A1 - A3) + DF_dG*(dA1_dx - dA3_dx);
   result_uoamput(3) = dDFdG_dF*(A1 - A3) + DF_dG*(dA1_dF - dA3_dF);
   result_uoamput(4) = d2DFdG_dx2*(A1 - A3) + 2*dDFdG_dx*(dA1_dx - dA3_dx) + DF_dG*(d2A1_dx2 - d2A3_dx2);
   result_uoamput(5) = d2DFdG_dF2*(A1 - A3) + 2*dDFdG_dF*(dA1_dF - dA3_dF) + DF_dG*(d2A1_dF2 - d2A3_dF2);
   result_uoamput(6) = dDFdG_dt*(A1 - A3) + DF_dG*(dA1_dt - dA3_dt);
   result_uoamput(7) = dDFdG_dsigma*(A1 - A3) + DF_dG*(dA1_dsigma - dA3_dsigma);
 endif
 if (result_uoamput(1)<=0)
    result_uoamput = zeros(1,7);
 endif 
endfunction

function [result_diamput] = diamput(F_bid,F_ask,barrier,strike,issue_date,expire_date,PPO,OSO,type)
 global DCC;
 global FC_DOM;
 global DSQ_Bid;
 global DSQ_Ask;
 global DSB_Bid;
 global DSB_Ask;
 tau = year_frac(issue_date,expire_date,DCC);
 sigma = ImpVol(strike,issue_date,expire_date);
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
   F = F_ask;
   DF_d = DF_d_ask;
   DF_f = DF_f_bid;
  else
   DF_dG = DF_d_bid_mod;
   F = F_bid;
   DF_d = DF_d_bid;
   DF_f = DF_f_ask;
 endif
 result_DFdG = differ_DFd(DF_dG,tau_mod);
 dDFdG_dx = result_DFdG(1);
 dDFdG_dF = result_DFdG(2);
 d2DFdG_dx2 = result_DFdG(3);
 d2DFdG_dF2 = result_DFdG(4);
 dDFdG_dt = result_DFdG(5);
 dDFdG_dsigma = result_DFdG(6);
 if (barrier < strike)
   result_A2 = differ_A2(F,DF_d,DF_f,tau,sigma,barrier,strike,-1);
   A2 = result_A2(1);
   dA2_dx = result_A2(2);
   dA2_dF = result_A2(3);
   d2A2_dx2 = result_A2(4);
   d2A2_dF2 = result_A2(5);
   dA2_dt = result_A2(6);
   dA2_dsigma = result_A2(7);
   result_A3 = differ_A3(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,1);
   A3 = result_A3(1);
   dA3_dx = result_A3(2);
   dA3_dF = result_A3(3);
   d2A3_dx2 = result_A3(4);
   d2A3_dF2 = result_A3(5);
   dA3_dt = result_A3(6);
   dA3_dsigma = result_A3(7);
   result_A4 = differ_A4(F,DF_d,DF_f,tau,sigma,barrier,strike,-1,1);
   A4 = result_A4(1);
   dA4_dx = result_A4(2);
   dA4_dF = result_A4(3);
   d2A4_dx2 = result_A4(4);
   d2A4_dF2 = result_A4(5);
   dA4_dt = result_A4(6);
   dA4_dsigma = result_A4(7);
   result_diamput(1) = DF_dG*(A2 - A3 + A4);
   result_diamput(2) = dDFdG_dx*(A2 - A3 + A4) + DF_dG*(dA2_dx - dA3_dx + dA4_dx);
   result_diamput(3) = dDFdG_dF*(A2 - A3 + A4) + DF_dG*(dA2_dF - dA3_dF + dA4_dF);
   result_diamput(4) = d2DFdG_dx2*(A2 - A3 + A4) + 2*dDFdG_dx*(dA2_dx - dA3_dx + dA4_dx) + DF_dG*(d2A2_dx2 - d2A3_dx2 + d2A4_dx2);
   result_diamput(5) = d2DFdG_dF2*(A2 - A3 + A4) + 2*dDFdG_dF*(dA2_dF - dA3_dF + dA4_dF) + DF_dG*(d2A2_dF2 - d2A3_dF2 + d2A4_dF2);
   result_diamput(6) = dDFdG_dt*(A2 - A3 + A4) + DF_dG*(dA2_dt - dA3_dt + dA4_dt);
   result_diamput(7) = dDFdG_dsigma*(A2 - A3 + A4) + DF_dG*(dA3_dsigma - dA3_dsigma + dA4_dsigma);
  else
   result_A1 = differ_A1(F,DF_d,DF_f,tau,sigma,strike,-1);
   A1 = result_A1(1);
   dA1_dx = result_A1(2);
   dA1_dF = result_A1(3);
   d2A1_dx2 = result_A1(4);
   d2A1_dF2 = result_A1(5);
   dA1_dt = result_A1(6);
   dA1_dsigma = result_A1(7);
   result_diamput(1) = DF_dG*A1;
   result_diamput(2) = dDFdG_dx*A1 + DF_dG*dA1_dx;
   result_diamput(3) = dDFdG_dF*A1 + DF_dG*dA1_dF;
   result_diamput(4) = d2DFdG_dx2*A1 + 2*dDFdG_dx*dA1_dx + DF_dG*d2A1_dx2;
   result_diamput(5) = d2DFdG_dF2*A1 + 2*dDFdG_dF*dA1_dF + DF_dG*d2A1_dF2;
   result_diamput(6) = dDFdG_dt*A1 + DF_dG*dA1_dt;
   result_diamput(7) = dDFdG_dsigma*A1 + DF_dG*dA1_dsigma;
 endif
 if (result_diamput(1)<=0)
    result_diamput = zeros(1,7);
 endif 
endfunction
