# domestic discount factor table creation
1;



%interpolacje do wyboru:
%linear on df
%linear on rate
%raw
%linear on lograte
%nat cubic
%Bspline
interp_method='linear on lograte';


 

%wygenerowanie tablic DS
[DSD_Bid, DSD_Ask, DSD_Ave]=DOM_curve_constr;

% zapis danych potrzebnych dla inicjalizacji rynku bez obliczen

%init_dom_data_gen();
 