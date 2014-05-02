% TO WERSJA DO PORTALU
% USUNIETE FUNKCJE KTORE NIE SA WYKORZYSTYWANE
% PRZY PRACY W PORTALU


%FINANSE OBLICZENIOWE
%DUŻY PROJEKT 2011
%FUNKCJE DYSKONTOWE
%ADAM RYTERSKI
%---------------------------
%POPRAWKI I MODYFIKACJE
%DUŻY PROJEKT 2012
%MARCIN SOSNOWSKI
%---------------------------
%POPRAWKI I MODYFIKACJE
%DUŻY PROJEKT 2013
%Bartosz Mielczarek
% obecna wersja:
%1. wprowadzona obsluga kwotowan futures - dane kwotowan futures zapisywane jako
%  FUT_DCC = FRA_DCC;%  FUT_BDA = FRA_BDA;%  FUT_EMA = FRA_EMA;%  FUT_rates = FRA_rates;
% uwaga w FUT_rates tenor zapisywane jako Mar10 (bez spacji)
%2. wprowadzone obcinanie za dlugich danych
% przyjeto zasade, ze ubcina sie kwotowania wczesniejszego instrumentu
% tzn. jesli Depo jest do 1Y a FRA zaczynaja sie od 3x6, to do obliczen
% bierze sie Depo do 3M a potem FRA. 
% analogicznie z FRA i IRS lub Futures i IRS
%3. zmodyfikowano uzycie interpolacji "nat cubic". Interpolacja jest na czynnikach dyskontowych
%4. wprowadzono interpolacje "Bspline"

1;
%##################GENEROWANIE market_ir/market_fx################################

%Funkcja 1------------------------------------------------------------------------
%funkcja generuje tablice market_ir
function [market_ir]=market_ir_gen()

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


% dane o kontraktach futures zapisywane sa pod zmiennymi odpowiadajacymi FRA
% dletego nalezy je przypisac zmiennym FUT
if (strcmp(CURR_FRAFUT,"Futures")==1) 
   FUT_DCC = FRA_DCC;   FUT_BDA = FRA_BDA;   FUT_EMA = FRA_EMA;   FUT_rates = FRA_rates; 
end   

n1=length(Depo_rates);
if (strcmp(CURR_FRAFUT,"FRA")==1)
    n2=length(FRA_rates);
else
    n2=length(FUT_rates);
endif
n3=length(IRS_rates);

 

%najpierw tworzymy tablice bez sprawdzania czy mamy ON/TN
market_ir2=cell(n1+n2+n3,7);

%DEP
for i=1:n1
  market_ir2{i,1}="DEP";
  market_ir2{i,2}=Depo_rates{i,1};
  [pay_dt,pay_fr]=cpn_dt(start_date,Depo_rates{i,1},"DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
  market_ir2{i,3}=pay_dt{1};
  market_ir2{i,4}=pay_dt{2}; 
  market_ir2{i,5}=Depo_rates{i,2}; %stopa BID
  market_ir2{i,6}=Depo_rates{i,3}; %stopa ASK
  market_ir2{i,7}=Depo_DCC; %konwencja
end

if (strcmp(CURR_FRAFUT,"FRA")==1) %FRA
for i=1:n2
  market_ir2{i+n1,1}="FRA";
  market_ir2{i+n1,2}=FRA_rates{i,1};
  [pay_dt,pay_fr]=cpn_dt(start_date,FRA_rates{i,1},"FRA",FC_DOM,FRA_DCC,FRA_BDA,FRA_day_to_spot,FRA_EMA);
  market_ir2{i+n1,3}=pay_dt{2};
  market_ir2{i+n1,4}=pay_dt{3}; 
  market_ir2{i+n1,5}=FRA_rates{i,2}; %stopa BID
  market_ir2{i+n1,6}=FRA_rates{i,3}; %stopa ASK
  market_ir2{i+n1,7}=FRA_DCC; %konwencja
end
else  %FUT
for i=1:n2
  market_ir2{i+n1,1}="FUT";
  market_ir2{i+n1,2}=FUT_rates{i,1};
  [pay_dt,pay_fr]=cpn_dt(start_date,FUT_rates{i,1},"FUT",FC_DOM,FUT_DCC,FUT_BDA,0,FUT_EMA);
  market_ir2{i+n1,3}=pay_dt{1};
  market_ir2{i+n1,4}=pay_dt{2}; 
  market_ir2{i+n1,5}=100-FUT_rates{i,2}-FUT_rates{i,3}; %stopa BID
  market_ir2{i+n1,6}=100-FUT_rates{i,2}-FUT_rates{i,3}; %stopa ASK
  market_ir2{i+n1,7}=FUT_DCC; %konwencja
end
end

%IRS
%daty tylko dla najdluzszego IRS
[pay_dt,pay_fr]=cpn_dt(start_date,IRS_rates{n3,1},"IRS",FC_DOM,IRS_DCC,IRS_BDA,IRS_day_to_spot,IRS_EMA,IRS_cf);
for i=1:n3
  tenor=IRS_rates{i,1};
  a=str2num(tenor(1:end-1)); %dlugosc IRS
  market_ir2{i+n1+n2,1}="IRS";
  market_ir2{i+n1+n2,2}=IRS_rates{i,1}; %tenor
  market_ir2{i+n1+n2,3}=pay_dt{1}; %początek
  market_ir2{i+n1+n2,4}=pay_dt{a*IRS_cf+1};%koniec
  market_ir2{i+n1+n2,5}=IRS_rates{i,2};%bid
  market_ir2{i+n1+n2,6}=IRS_rates{i,3};%ask
  market_ir2{i+n1+n2,7}=IRS_DCC;%konwencja
end

% tu dokonuje sie obcinania zachodzacych danych
% obcinanie danych dotyczy zawsze najdluzszych zapadalnosci danego instrumentu

indic = ones(1,n1+n2+n3);

for i=2:(n1+n2+n3)
    
  if day_diff(market_ir2{i-1,4}, market_ir2{i,4}, "ACT") <= 0
     for k=1:(i-1)
        if day_diff(market_ir2{i-k,4}, market_ir2{i,4}, "ACT") <= 0  
           indic(i-k) = 0;
        end   
     end
  end
end

ni = sum(indic);

market_ir3=cell(ni,7);

k=1;
for i=1:(n1+n2+n3)
  if (indic(i) == 1)
     market_ir3{k,1}=market_ir2{i,1};
	  market_ir3{k,2}=market_ir2{i,2};
	  market_ir3{k,3}=market_ir2{i,3};
	  market_ir3{k,4}=market_ir2{i,4};
	  market_ir3{k,5}=market_ir2{i,5};
	  market_ir3{k,6}=market_ir2{i,6};
	  market_ir3{k,7}=market_ir2{i,7};
     k = k+1;
  end
end
     
     
% w tej funkcji zapisuje sie do market_ir_long   to co podal inwestor
% linie 174-193 mozna wykasowac !!!
%market_ir_long=market_ir2;

%zapis do pliku
%n=length(market_ir_long);
%file_id = fopen('market_irl.m', 'w');

%fprintf(file_id, "market_ir={");


%for i=1:n-1
%    fprintf(file_id, "\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, \"%s\";",market_ir_long{i,1},market_ir_long{i,2},...
%    market_ir_long{i,3},market_ir_long{i,4},market_ir_long{i,5},market_ir_long{i,6},market_ir_long{i,7});
%end

%fprintf(file_id, "\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, \"%s\"};",market_ir_long{n,1},market_ir_long{n,2},...
%    market_ir_long{n,3},market_ir_long{n,4},market_ir_long{n,5},market_ir_long{n,6},market_ir_long{n,7});

%fclose(file_id);


% w tej funkcji zapisuje sie do market_ir tylko to co podal inwestor z obcieciem zachodzacych danych

market_ir=market_ir3;

%zapis do pliku
n=length(market_ir);

file_id = fopen('market_ir.m', 'w');

fprintf(file_id, "market_ir={");

for i=1:n-1
    fprintf(file_id, "\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, \"%s\";",market_ir{i,1},market_ir{i,2},...
    market_ir{i,3},market_ir{i,4},market_ir{i,5},market_ir{i,6},market_ir{i,7});
end


fprintf(file_id, "\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, \"%s\"};",market_ir{n,1},market_ir{n,2},...
market_ir{n,3},market_ir{n,4},market_ir{n,5},market_ir{n,6},market_ir{n,7});

fclose(file_id);


endfunction


%Funkcja 2------------------------------------------------------------------------
%funkcja generuje tablice market_fx
function [market_fx]=market_fx_gen()
global FC_FOR;
global start_date;

global FOR_DCC;
global FOR_BDA;
global FOR_EMA;
global FOR_day_to_spot;
global FX_rate;
global SWAP_POINTS;
global PDR;
global pip_val;
global CURR_DOM;
global CURR_FOR;
global BAS_CURR;
global QUO_CURR;


n=length(SWAP_POINTS);

market_fx=cell(n+1,4);

%sztucznie do tablicy dodajemy pierwszy wiersz z datą start_date oraz kursami spot
%Taki zabieg eliminuje problem ekstrapolacji w "lewą stronę"
market_fx{1,1}=start_date; 
market_fx{1,2}=FX_rate(1);
market_fx{1,3}=FX_rate(2);
market_fx{1,4}=FOR_DCC;

BAS_CURR =  base_curr(CURR_FOR, CURR_DOM);

if BAS_CURR == CURR_DOM
   QUO_CURR = CURR_FOR;
else
    QUO_CURR = CURR_DOM;
end       

for i=1:n
  [pay_dt,pay_fr]=cpn_dt(start_date,SWAP_POINTS{i,1},"FWD",FC_FOR,FOR_DCC,FOR_BDA,FOR_day_to_spot,FOR_EMA);
  market_fx{i+1,1}=pay_dt{2};%data
  market_fx{i+1,2}=PDR*abs(SWAP_POINTS{i,2})*pip_val+FX_rate(1); %stopa BID
  market_fx{i+1,3}=PDR*abs(SWAP_POINTS{i,3})*pip_val+FX_rate(2); %stopa ASK
  market_fx{i+1,4}=FOR_DCC; %konwencja
end




%zapis do pliku
n=length(market_fx);
file_id = fopen('market_fxr.m', 'w');

 fprintf(file_id, "BAS_CURR = \"%s\";\n",  BAS_CURR);
 fprintf(file_id, "QUO_CURR = \"%s\";\n\n",  QUO_CURR);

fprintf(file_id, "market_fx={");


for i=1:n-1
    fprintf(file_id, "\"%s\", %f, %f, \"%s\";",market_fx{i,1},market_fx{i,2},...
    market_fx{i,3},market_fx{i,4});
end


fprintf(file_id, "\"%s\", %f, %f, \"%s\"};",market_fx{n,1},market_fx{n,2},...
market_fx{n,3},market_fx{n,4});

fclose(file_id);

endfunction


%Funkcja 3------------------------------------------------------------------------
%funkcja zwraca daty przeplywow pien. dla konkretnych instrumentow
function [pay_dt,pay_fr]=cpn_dt(start_date,one_tenor,market,FC,DCC,BDA,Day_to_spot,EMA,cf)

if nargin<8 %gdy nie potrzebujemy cf = liczba kuponow IRS w roku
  cf=1;
end

EMA=str2num(EMA);

switch (market)
    case "FRA"
	n=index(one_tenor,"X"); %pozycja X
	a=str2num(one_tenor(1:n-1)); %ilosc miesiecy do poczatku kontraktu
	b=str2num(one_tenor(n+1:length(one_tenor))); %ilosc miesiecy do konca kontraktu
	[pay_dt,pay_fr]=payments4(start_date,[a,b-a],FC,DCC,BDA,Day_to_spot,EMA);

    case "FUT"
   if str2num(substr(one_tenor,4,2))<40 
	   a=["01-",substr(one_tenor,1,3),"-20",substr(one_tenor,4,2)]; % jesli daty od roku 2000 i one_tenor = "Mar10" to a = "01-Mar-2010"
	else 
	   a=["01-",substr(one_tenor,1,3),"-19",substr(one_tenor,4,2)]; % jesli daty przed rokiem 2000
	end   
	[m1,m2]=wedn(a,one_tenor); % znajduje 3 Wednesday w miesiacu one_tenor
	b=[m1,"-",m2,"-",substr(a,8,4)]; % data 3 Wednesday
	c=year_frac(start_date,b,DCC)/12; %integer (ilosc dni do poczatku kontraktu od start_date)
	[pay_dt,pay_fr]=payments4(b,3,FC,DCC,BDA,Day_to_spot,EMA);
	
    case "IRS"
	n=index(one_tenor,"Y"); %pozycja Y
	a=str2num(one_tenor(1:n-1)); %ilosc lat
	[pay_dt,pay_fr]=payments5(start_date,cf*a,12/cf,FC,DCC,BDA,Day_to_spot,EMA);	

    case {"DEP", "FWD"}
	  T=one_tenor;
	      	switch (T)
		    case "ON"
		      [pay_dt,pay_fr]=payments(start_date,1,FC,DCC,"sfbd",0);
		    case "TN"
		      [pay_dt,pay_fr]=payments(start_date,1,FC,DCC,"sfbd",1);
		    case "SN"
		      [pay_dt,pay_fr]=payments(start_date,1,FC,DCC,BDA,Day_to_spot);
		    case {"SW","1W"}
		      [pay_dt,pay_fr]=payments(start_date,7,FC,DCC,BDA,Day_to_spot);
		    case "2W"
		      [pay_dt,pay_fr]=payments(start_date,14,FC,DCC,BDA,Day_to_spot);
		    case "3W"
		      [pay_dt,pay_fr]=payments(start_date,21,FC,DCC,BDA,Day_to_spot);
		    case "4W"
		      [pay_dt,pay_fr]=payments4(start_date,1,FC,DCC,BDA,Day_to_spot,EMA);
		    otherwise
		      if index(T,"M")>0 %jeśli mamy miesiące
			[pay_dt,pay_fr]=payments4(start_date,str2num(T(1:end-1)),FC,DCC,BDA,Day_to_spot,EMA);
		      elseif index(T,"Y")>0 %jeśli lata
			shiftM=12*str2num(T(1:end-1)); %lata zamieniamy na miesiące
			[pay_dt,pay_fr]=payments4(start_date,shiftM,FC,DCC,BDA,Day_to_spot,EMA);
		      end
		end %switch
end %switch

endfunction

%#####################DS dla waluty krajowej######################################

%Funkcje 4 i 5--------------------------------------------------------------------
%funkcje pomocnicze
%zamieniajace stope proc. na df i odwrotnie
%t - ulamek roku
function [df]=rate2df(t,rate, method)
t=t(:);
rate=rate(:);

if nargin<3 %gdy nie mamy metody
  method="continuous";
end

switch (method)
  case "continuous"
    df = exp(-t.*rate);
  case "simple"
    df = (t.*rate+1).^(-1);
end

df=df(:);

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rate]=df2rate(t,df, method)
t=t(:);
df=df(:);

if nargin<3 %gdy nie mamy metody
  method="continuous";
end


switch (method)
  case "continuous"
    rate = -log(df)./t;
  case "simple"
    rate = (df.^(-1)-1)./t;
end

rate=rate(:);

endfunction


%Funkcja 6------------------------------------------------------------------------
%funkcja sluzaca do budowy tablic czynnikow dyskontowych dla waluty krajowej
function [DSD_Bid, DSD_Ask, DSD_Ave]=DOM_curve_constr()

global IRS_cf;
global IRS_BDA;
global IRS_EMA;
global FC_DOM;

global Depo_DCC;
global Depo_BDA;
global Depo_day_to_spot;
global Depo_EMA;
global Depo_rates;

global start_date;
global interp_method;

% For discounting during the process of discount table construction
% we cannot use spline interpolation
% Hence if interp_method is one of spline interpolations
% we replece it by raw_interp
% Below we assign to logical variable simple_interp value TRUE if interp_method is not spline interpolation

simple_int_m = {"linear on df", "linear on rate", "raw", "linear on lograte"};
wsk = 0;
for i = 1 : columns(simple_int_m)
if strcmp(simple_int_m(i), interp_method) == 1
wsk = i;
end  
end

simple_interp = false;

if wsk > 0
simple_interp = true;
end




source market_ir.m;

market_ir2 = market_ir;

n=length(market_ir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sprawdzenie czy mamy ON i TN
if Depo_day_to_spot==0 %nie potrzebujemy ON
    market_ir=market_ir2;

elseif Depo_day_to_spot==1 %TN=SN
    if strcmp(Depo_rates{1,1},"ON")==1 %gdy mamy ON
      market_ir=market_ir2;


    else %gdy nie mamy ON
      market_ir=cell(n+1,7); %dodajemy dodatkowy wiersz do market_ir
      %------ON----------
      market_ir{1,1}="DEP";
      market_ir{1,2}="ON";
      %daty
      [pay_dt,pay_fr]=cpn_dt(start_date,"ON","DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
      market_ir{1,3}=pay_dt{1};
      market_ir{1,4}=pay_dt{2};
      %nie mamy jeszcze wsp. wezlow
      x(1)=year_frac(start_date,market_ir2{1,4},Depo_DCC);
      x(2)=year_frac(start_date,market_ir2{2,4},Depo_DCC);
      %interpolacja ON
      market_ir{1,5}=linear_interp(pay_fr{1},x,[market_ir2{1,5},market_ir2{2,5}]); %stopa BID
      market_ir{1,6}=linear_interp(pay_fr{1},x,[market_ir2{1,6},market_ir2{2,6}]); %stopa ASK
      %konwencja
      market_ir{1,7}=Depo_DCC; 
      
      
      for i=1:length(market_ir2)
	  market_ir{i+1,1}=market_ir2{i,1};
	  market_ir{i+1,2}=market_ir2{i,2};
	  market_ir{i+1,3}=market_ir2{i,3};
	  market_ir{i+1,4}=market_ir2{i,4};
	  market_ir{i+1,5}=market_ir2{i,5};
	  market_ir{i+1,6}=market_ir2{i,6};
	  market_ir{i+1,7}=market_ir2{i,7};
      end
    end


elseif Depo_day_to_spot>=2 %gdy potrzebujemy ON i TN
    if strcmp(Depo_rates{1,1},"ON")==1 && strcmp(Depo_rates{2,1},"TN")==1 %mamy ON i TN
      market_ir=market_ir2;

    elseif strcmp(Depo_rates{1,1},"ON")==1 && strcmp(Depo_rates{2,1},"TN")==0 %mamy ON, nie mamy TN
      market_ir=cell(n+1,7); %dodajemy dodatkowy wiersz do market_ir

      %----ON----
      market_ir{1,1}=market_ir2{1,1};
      market_ir{1,2}=market_ir2{1,2};
      market_ir{1,3}=market_ir2{1,3};
      market_ir{1,4}=market_ir2{1,4};
      market_ir{1,5}=market_ir2{1,5};
      market_ir{1,6}=market_ir2{1,6};
      market_ir{1,7}=market_ir2{1,7};  
      %----TN----
      market_ir{2,1}="DEP";
      market_ir{2,2}="TN";
      %daty
      [pay_dt,pay_fr]=cpn_dt(start_date,"TN","DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
      market_ir{2,3}=pay_dt{1};
      market_ir{2,4}=pay_dt{2};
      %nie mamy jeszcze wsp. wezlow
      x(1)=year_frac(start_date,market_ir2{1,4},Depo_DCC);
      x(2)=year_frac(start_date,market_ir2{2,4},Depo_DCC);
      %interpolacja TN
      market_ir{2,5}=linear_interp(x(1)+pay_fr{1},x,[market_ir2{1,5},market_ir2{2,5}]); %stopa BID
      market_ir{2,6}=linear_interp(x(1)+pay_fr{1},x,[market_ir2{1,6},market_ir2{2,6}]); %stopa ASK
      %konwencja
      market_ir{2,7}=Depo_DCC; 

      for i=2:length(market_ir2)
	  market_ir{i+1,1}=market_ir2{i,1};
	  market_ir{i+1,2}=market_ir2{i,2};
	  market_ir{i+1,3}=market_ir2{i,3};
	  market_ir{i+1,4}=market_ir2{i,4};
	  market_ir{i+1,5}=market_ir2{i,5};
	  market_ir{i+1,6}=market_ir2{i,6};
	  market_ir{i+1,7}=market_ir2{i,7};
      end

    elseif strcmp(Depo_rates{1,1},"ON")==0 && strcmp(Depo_rates{1,1},"TN")==1 %mamy TN, nie mamy ON
      market_ir=cell(n+1,7); %dodajemy dodatkowy wiersz do market_ir

      %----TN----
      market_ir{2,1}=market_ir2{1,1};
      market_ir{2,2}=market_ir2{1,2};
      market_ir{2,3}=market_ir2{1,3};
      market_ir{2,4}=market_ir2{1,4};
      market_ir{2,5}=market_ir2{1,5};
      market_ir{2,6}=market_ir2{1,6};
      market_ir{2,7}=market_ir2{1,7};  
      %----ON----
      market_ir{1,1}="DEP";
      market_ir{1,2}="ON";
      %daty
      [pay_dt,pay_fr]=cpn_dt(start_date,"ON","DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
      market_ir{1,3}=pay_dt{1};
      market_ir{1,4}=pay_dt{2};
      %nie mamy jeszcze wsp. wezlow
      x(1)=year_frac(start_date,market_ir2{1,4},Depo_DCC); %koniec TN
      x(2)=year_frac(start_date,market_ir2{2,4},Depo_DCC);  %koniec nastepnego DEP
      %interpolacja ON
      market_ir{1,5}=linear_interp(pay_fr{1},x,[market_ir2{1,5},market_ir2{2,5}]); %stopa BID
      market_ir{1,6}=linear_interp(pay_fr{1},x,[market_ir2{1,6},market_ir2{2,6}]); %stopa ASK
      %konwencja
      market_ir{1,7}=Depo_DCC; 

      for i=2:length(market_ir2)
	  market_ir{i+1,1}=market_ir2{i,1};
	  market_ir{i+1,2}=market_ir2{i,2};
	  market_ir{i+1,3}=market_ir2{i,3};
	  market_ir{i+1,4}=market_ir2{i,4};
	  market_ir{i+1,5}=market_ir2{i,5};
	  market_ir{i+1,6}=market_ir2{i,6};
	  market_ir{i+1,7}=market_ir2{i,7};
      end

    elseif strcmp(Depo_rates{1,1},"ON")==0 && strcmp(Depo_rates{2,1},"TN")==0 %nie mamy ON, nie mamy TN
      market_ir=cell(n+2,7); %dodajemy 2 dodatkowe wiersze do market_ir

      %----ON----
      market_ir{1,1}="DEP";
      market_ir{1,2}="ON";
      %daty
      [pay_dt,pay_frO]=cpn_dt(start_date,"ON","DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
      market_ir{1,3}=pay_dt{1};
      market_ir{1,4}=pay_dt{2};
      %nie mamy jeszcze wsp. wezlow
      x(1)=year_frac(start_date,market_ir2{1,4},Depo_DCC);
      x(2)=year_frac(start_date,market_ir2{2,4},Depo_DCC);
      %interpolacja ON
      market_ir{1,5}=linear_interp(pay_frO{1},x,[market_ir2{1,5},market_ir2{2,5}]); %stopa BID
      market_ir{1,6}=linear_interp(pay_frO{1},x,[market_ir2{1,6},market_ir2{2,6}]); %stopa ASK
      %konwencja
      market_ir{1,7}=Depo_DCC; 

      %----TN----
      market_ir{2,1}="DEP";
      market_ir{2,2}="TN";
      %daty
      [pay_dt,pay_frT]=cpn_dt(start_date,"TN","DEP",FC_DOM,Depo_DCC,Depo_BDA,Depo_day_to_spot,Depo_EMA);
      market_ir{2,3}=pay_dt{1};
      market_ir{2,4}=pay_dt{2};
      %nie mamy jeszcze wsp. wezlow
      x(1)=year_frac(start_date,market_ir2{1,4},Depo_DCC);
      x(2)=year_frac(start_date,market_ir2{2,4},Depo_DCC);
      %interpolacja TN
      market_ir{2,5}=linear_interp(pay_frO{1}+pay_frT{1},x,[market_ir2{1,5},market_ir2{2,5}]); %stopa BID
      market_ir{2,6}=linear_interp(pay_frO{1}+pay_frT{1},x,[market_ir2{1,6},market_ir2{2,6}]); %stopa ASK
      %konwencja
      market_ir{2,7}=Depo_DCC; 


      for i=1:length(market_ir2)
	  market_ir{i+2,1}=market_ir2{i,1};
	  market_ir{i+2,2}=market_ir2{i,2};
	  market_ir{i+2,3}=market_ir2{i,3};
	  market_ir{i+2,4}=market_ir2{i,4};
	  market_ir{i+2,5}=market_ir2{i,5};
	  market_ir{i+2,6}=market_ir2{i,6};
	  market_ir{i+2,7}=market_ir2{i,7};
      end



    end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(market_ir);
zzz=0;

DSD_Bid=cell(n+1,3);
DSD_Ask=cell(n+1,3);

%dodajemy sztucznie pierwszy wiersz, aby uniknac problemu z ekstrapolacja w "lewe strone"
DSD_Bid(1,1)=start_date;
DSD_Bid(1,2)=1;
DSD_Bid(1,3)="ACT/365";
DSD_Ask(1,1)=start_date;
DSD_Ask(1,2)=1;
DSD_Ask(1,3)="ACT/365";
%tablice pomocnicze zawierające czynnik dyskontowy i ulamek roku, jw dodajemy sztuczny pierwszy wiersz
tDSA(1,1)=0;
tDSA(1,2)=1;
tDSB(1,1)=0;
tDSB(1,2)=1;


%w zaleznosci od depo_day_to_spot liczymy czynnik dyskontowy spot
if Depo_day_to_spot==0
    spot=0;
    df_spotA=1;
    df_spotB=1;
    st=1; %zmienna pomocnicza, okreslajaca miejsce startu glownej petli

elseif Depo_day_to_spot==1
    spot=year_frac(start_date,market_ir{1,4},market_ir{1,7});%koniec ON to spot
    df_spotB=rate2df(spot,market_ir{1,5}/100,"simple");
    df_spotA=rate2df(spot,market_ir{1,6}/100,"simple");

    tDSA(2,1)=spot;
    tDSA(2,2)=df_spotA;
    tDSB(2,1)=spot;
    tDSB(2,2)=df_spotB;    

    DSD_Bid(2,1)=market_ir{1,4};
    DSD_Bid(2,2)=df_spotB;
    DSD_Bid(2,3)=market_ir{1,7};
    DSD_Ask(2,1)=market_ir{1,4};
    DSD_Ask(2,2)=df_spotA;
    DSD_Ask(2,3)=market_ir{1,7};
    st=2;

elseif Depo_day_to_spot>=2
    %-----ON-----
    tDSB(2,1)=year_frac(start_date,market_ir{1,4},market_ir{1,7});%koniec ON
    tDSB(2,2)=rate2df(tDSB(2,1),market_ir{1,5}/100,"simple");
    tDSA(2,1)=tDSB(2,1);%koniec ON
    tDSA(2,2)=rate2df(tDSA(2,1),market_ir{1,6}/100,"simple");

    DSD_Bid(2,1)=market_ir{1,4};
    DSD_Bid(2,2)=tDSB(2,2);
    DSD_Bid(2,3)=market_ir{1,7};
    DSD_Ask(2,1)=market_ir{1,4};
    DSD_Ask(2,2)=tDSA(2,2);
    DSD_Ask(2,3)=market_ir{1,7};

    %-----TN-----
    tDSB(3,1)=year_frac(start_date,market_ir{2,4},market_ir{2,7});%koniec TN
    tDSB(3,2)=rate2df(tDSB(3,1)-tDSB(2,1),market_ir{2,5}/100,"simple")*tDSB(2,2);
    tDSA(3,1)=tDSB(3,1);%koniec TN
    tDSA(3,2)=rate2df(tDSA(3,1)-tDSA(2,1),market_ir{2,6}/100,"simple")*tDSA(2,2);

    DSD_Bid(3,1)=market_ir{2,4};
    DSD_Bid(3,2)=tDSB(3,2);
    DSD_Bid(3,3)=market_ir{2,7};
    DSD_Ask(3,1)=market_ir{2,4};
    DSD_Ask(3,2)=tDSA(3,2);
    DSD_Ask(3,3)=market_ir{2,7};

    %liczymy czynnik dyskontowy spot
    spot=year_frac(start_date,market_ir{3,3},market_ir{3,7});%3 depozyt zaczyna sie od spot
    
    if simple_interp
       df_spotB=addinterp(spot,tDSB);%df_spot liczymy z ekstrapolacji, mamy 3 wezly wiec ok
       df_spotA=addinterp(spot,tDSA); 
    else
       df_spotB=raw_interp(spot,tDSB(:,1),tDSB(:,2)); % tak liczymy gdy intepolacja splinowa
       df_spotA=raw_interp(spot,tDSA(:,1),tDSA(:,2)); 
    end   

    st=3;
end

%liczymy pozostale czynniki dyskontowe
for i=st:n
  switch (market_ir{i,1})
    case "DEP"
    %----------------BID----------------
	tDSB(i+1,1)=year_frac(start_date,market_ir{i,4},market_ir{i,7}); %odleglosc od start_date do konca dep.
	DF=df_spotB/((tDSB(i+1,1)-spot)*(market_ir{i,5})/100+1);%df
	tDSB(i+1,2)=DF;
	DSD_Bid{i+1,1}=market_ir{i,4};
	DSD_Bid{i+1,2}=tDSB(i+1,2);
	DSD_Bid{i+1,3}=market_ir{i,7};
    %----------------ASK----------------
	tDSA(i+1,1)=tDSB(i+1,1); %odleglosc od start_date do konca dep.
	DF=df_spotA/((tDSA(i+1,1)-spot)*(market_ir{i,6})/100+1);%df
	tDSA(i+1,2)=DF;
	DSD_Ask{i+1,1}=market_ir{i,4};
	DSD_Ask{i+1,2}=tDSA(i+1,2);
	DSD_Ask{i+1,3}=market_ir{i,7};
	
	
    case {"FRA", "FUT"}
    %----------------BID----------------
	temp=year_frac(start_date,market_ir{i,3},market_ir{i,7});
   if simple_interp
       DF1=addinterp(temp,tDSB);%interpolujemy DF1
   else
        DF1=raw_interp(temp,tDSB(:,1),tDSB(:,2)); % tak liczymy gdy intepolacja splinowa
   end   
   
	tDSB(i+1,1)=year_frac(start_date,market_ir{i,4},market_ir{i,7});
	DF2=DF1/((tDSB(i+1,1)-temp)*(market_ir{i,5})/100+1);%liczymy czynnik dyskontujący dla daty końca kontraktu
	tDSB(i+1,2)=DF2;
	DSD_Bid{i+1,1}=market_ir{i,4};
	DSD_Bid{i+1,2}=tDSB(i+1,2);
	DSD_Bid{i+1,3}=market_ir{i,7};
    %----------------ASK----------------
	temp=year_frac(start_date,market_ir{i,3},market_ir{i,7});
	if simple_interp
       DF1=addinterp(temp,tDSA);%interpolujemy DF1
   else
        DF1=raw_interp(temp,tDSA(:,1),tDSA(:,2)); % tak liczymy gdy intepolacja splinowa
   end 
        
	tDSA(i+1,1)=year_frac(start_date,market_ir{i,4},market_ir{i,7});
	DF2=DF1/((tDSA(i+1,1)-temp)*(market_ir{i,6})/100+1);%liczymy czynnik dyskontujący dla daty końca kontraktu
	tDSA(i+1,2)=DF2;
	DSD_Ask{i+1,1}=market_ir{i,4};
	DSD_Ask{i+1,2}=tDSA(i+1,2);
	DSD_Ask{i+1,3}=market_ir{i,7};
	
	
	case "IRS"
	%-----------------------------------
	zzz=zzz+1;%zliczam kontrakty IRS
	tenor=market_ir{i,2};
	dlugosc=index(tenor,"Y");
	a=str2num(tenor(1:dlugosc-1));
	
	%konstrukcja tabeli pomocniczej zawierajacej informacje na temat kontraktow IRS
	IRS{zzz,1}=a;
	IRS{zzz,2}=market_ir{i,3};
	IRS{zzz,3}=market_ir{i,4};
	IRS{zzz,4}=0.01*market_ir{i,5};
	IRS{zzz,5}=0.01*market_ir{i,6};
	IRS{zzz,6}=market_ir{i,7};
	
   end %switch
end %for
	
	
[cashflows_dt,cashflows_fr]=payments5(IRS{zzz,2},IRS_cf*IRS{zzz,1},12/IRS_cf,FC_DOM,market_ir{i,7},IRS_BDA,0,str2num(IRS_EMA));%przeplywy aktywow kontraktow IRS
cashflows_totalfr=spot+cumsum(cell2mat(cashflows_fr));%daty przeplywow zeskalowane do lat
%1. kolumna - czas w latach, 2. kurs Bid, 3. kurs Ask
IRS_maturityfr=[];
for i=1:zzz
	IRS_maturityfr=[IRS_maturityfr;spot+year_frac(IRS{i,2},IRS{i,3},IRS{i,6})];
end %for
 
obrot=0;
%inicjalizacja algorytmu
DFstrzal=[[tDSB';tDSA(:,2)'],[IRS_maturityfr';[rate2df(IRS_maturityfr,[IRS{:,4}],"simple"),rate2df(IRS_maturityfr,[IRS{:,5}],"simple")]']];
DFinterp=[cashflows_totalfr,addinterp(cashflows_totalfr,DFstrzal(1:2,:)'),addinterp(cashflows_totalfr,[DFstrzal(1,:)',DFstrzal(3,:)'])];
Rx=[];
for i=1:zzz
	czastrwania=IRS{i,1};
	Rx=[Rx,[czastrwania;(1-DFinterp(czastrwania*IRS_cf,2))/(sum(cell2mat(cashflows_fr(1:czastrwania*IRS_cf)).*DFinterp(1:czastrwania*IRS_cf,2)));(1-DFinterp(czastrwania*IRS_cf,3))/(sum(cell2mat(cashflows_fr(1:czastrwania*IRS_cf)).*DFinterp(1:czastrwania*IRS_cf,3)))]];
end

%algorytm poszukiwania krzywej dyskontowej z kontraktow IRS, Hagan&West s.92
while (max(max(abs(Rx(2,:)'-cell2mat(IRS(:,4)))),max(abs(Rx(3,:)'-cell2mat(IRS(:,5))))) > eps) && (obrot<30) %ograniczenie zwiazane ze zbyt duza dokladnoscia
	DFstrzal=[tDSB';tDSA(:,2)'];%nowa propozycja czynnikow dyskontowych
	for i=1:zzz
		czastrwania=IRS{i,1};
		DFstrzal=[DFstrzal,[IRS_maturityfr(i);(1-IRS{i,4}*sum(cell2mat(cashflows_fr(1:(czastrwania*IRS_cf-1))).*DFinterp(1:(czastrwania*IRS_cf-1),2)))/(1+IRS{i,4}*cashflows_fr{czastrwania*IRS_cf});(1-IRS{i,5}*sum(cell2mat(cashflows_fr(1:(czastrwania*IRS_cf-1))).*DFinterp(1:(czastrwania*IRS_cf-1),3)))/(1+IRS{i,5}*cashflows_fr{czastrwania*IRS_cf})]];
	end %for
	DFinterp=[cashflows_totalfr,addinterp(cashflows_totalfr,DFstrzal(1:2,:)'),addinterp(cashflows_totalfr,[DFstrzal(1,:)',DFstrzal(3,:)'])];
	Rx=[];
	for i=1:zzz
		czastrwania=IRS{i,1};
		Rx=[Rx,[czastrwania;(1-DFinterp(czastrwania*IRS_cf,2))/(sum(cell2mat(cashflows_fr(1:czastrwania*IRS_cf)).*DFinterp(1:czastrwania*IRS_cf,2)));(1-DFinterp(czastrwania*IRS_cf,3))/(sum(cell2mat(cashflows_fr(1:czastrwania*IRS_cf)).*DFinterp(1:czastrwania*IRS_cf,3)))]];
	end %for
	obrot=obrot+1;
end %while
tDSB=DFstrzal(1:2,:)';
tDSA=[DFstrzal(1,:);DFstrzal(3,:)]';


%-----------------------------BID-------------------
for i=1:zzz
	DSD_Bid{n+1-zzz+i,1}=IRS{i,3};
	DSD_Bid{n+1-zzz+i,2}=tDSB(i+n+1-zzz,2);
	DSD_Bid{n+1-zzz+i,3}=IRS{i,6};
end %for

%------------------------------ASK------------------
for i=1:zzz
	DSD_Ask{n+1-zzz+i,1}=IRS{i,3};
	DSD_Ask{n+1-zzz+i,2}=tDSA(i+n+1-zzz,2);
	DSD_Ask{n+1-zzz+i,3}=IRS{i,6};
end %for

%czynniki dyskontowe sa juz obliczone

%zapis do pliku
n=length(DSD_Bid);
DSD_Ave = DSD_Ask;
for i=1:n
    DSD_Ave{i,2} = (DSD_Ask{i,2} + DSD_Bid{i,2})/2;
end

file_id = fopen('DSD_Bid.m', 'w');
file_id2 = fopen('DSD_Ask.m', 'w');
file_id3 = fopen('DSD_Ave.m', 'w');

fprintf(file_id, "DSD_Bid={");
fprintf(file_id2, "DSD_Ask={");
fprintf(file_id3, "DSD_Ave={");

for i=1:n-1
    fprintf(file_id, "\"%s\", %f, \"%s\";",DSD_Bid{i,1},DSD_Bid{i,2},DSD_Bid{i,3});
    fprintf(file_id2, "\"%s\", %f, \"%s\";",DSD_Ask{i,1},DSD_Ask{i,2},DSD_Ask{i,3});
    fprintf(file_id3, "\"%s\", %f, \"%s\";",DSD_Ave{i,1},DSD_Ave{i,2},DSD_Ave{i,3});
end


fprintf(file_id, "\"%s\", %f, \"%s\"};\n",DSD_Bid{n,1},DSD_Bid{n,2},DSD_Bid{n,3});
fprintf(file_id2, "\"%s\", %f, \"%s\"};\n",DSD_Ask{n,1},DSD_Ask{n,2},DSD_Ask{n,3});
fprintf(file_id3, "\"%s\", %f, \"%s\"};\n",DSD_Ave{n,1},DSD_Ave{n,2},DSD_Ave{n,3});
fclose(file_id);
fclose(file_id2);
fclose(file_id3);

endfunction

%################DS dla waluty zagranicznej#######################################

%Funkcja 7------------------------------------------------------------------------
%funkcja sluzaca do budowy tablic czynnikow dyskontowych dla waluty zagranicznej
function [DSF_Bid, DSF_Ask, DSF_Ave]=FOR_curve_constr()

global FX_rate;
global start_date;

global DSD_Bid;
global DSD_Ask;
global DSD_Ave;

global DSF_Ask;
global DSF_Bid;
global DSF_Ave;

global DSQ_Ask;
global DSQ_Bid;
global DSQ_Ave;

global DSB_Ask;
global DSB_Bid;
global DSB_Ave;

global CURR_DOM;
global CURR_FOR;
global BAS_CURR;
global QUO_CURR;
global CURR_PAIR;

source market_fxr.m;


n=length(market_fx);

DSF_Bid=cell(n,3);
DSF_Ask=cell(n,3);

ind =1;
     


if strcmp(start_date, market_fx(1,1))
    DSF_Bid{1,1}=market_fx{1,1};
    DSF_Ask{1,1}=market_fx{1,1};
    DSF_Bid{1,2}=1;
    DSF_Ask{1,2}=1;
    DSF_Bid{1,3}=market_fx{1,4};
    DSF_Ask{1,3}=market_fx{1,4};
    ind = 2;
endif     

disAsk = ones(1,n);
disBid = ones(1,n);

%liczymy czynniki dyskontujace dla dat z kwotowan punktow swapowych
for i = ind:n
	disAsk(i)=DF(start_date, market_fx(i,1),DSD_Ask);
	disBid(i)=DF(start_date, market_fx(i,1),DSD_Bid);
endfor

if (BAS_CURR==CURR_DOM & CURR_PAIR=="1") % jesli kurs wprowadzono BAS/QUO a BAS = DOM
  for i=ind:n
    DSF_Bid{i,1}=market_fx{i,1};
    DSF_Ask{i,1}=market_fx{i,1};
    DSF_Bid{i,2}=FX_rate(1)*disAsk(i)/market_fx{i,2};
    DSF_Ask{i,2}=FX_rate(2)*disBid(i)/market_fx{i,3};
    DSF_Bid{i,3}=market_fx{i,4};
    DSF_Ask{i,3}=market_fx{i,4};
  end
else % jesli kurs wprowadzono FOR/DOM 
  for i=ind:n
    DSF_Bid{i,1}=market_fx{i,1};
    DSF_Ask{i,1}=market_fx{i,1};
    DSF_Ask{i,2}=disBid(i)*market_fx{i,2}/FX_rate(1);
    DSF_Bid{i,2}=disAsk(i)*market_fx{i,3}/FX_rate(2);
    DSF_Bid{i,3}=market_fx{i,4};
    DSF_Ask{i,3}=market_fx{i,4};
  end
end 


n=length(DSF_Bid);
DSF_Ave = DSF_Ask;
for i=1:n
    DSF_Ave{i,2} = (DSF_Ask{i,2} + DSF_Bid{i,2})/2;
end

%zapis do pliku

file_id = fopen('DSF_Bid.m', 'w');
file_id2 = fopen('DSF_Ask.m', 'w');
file_id3 = fopen('DSF_Ave.m', 'w');

fprintf(file_id, "DSF_Bid={");
fprintf(file_id2, "DSF_Ask={");
fprintf(file_id3, "DSF_Ave={");

for i=1:n-1
    fprintf(file_id, "\"%s\", %f, \"%s\";",DSF_Bid{i,1},DSF_Bid{i,2},DSF_Bid{i,3});
    fprintf(file_id2, "\"%s\", %f, \"%s\";",DSF_Ask{i,1},DSF_Ask{i,2},DSF_Ask{i,3});
    fprintf(file_id3, "\"%s\", %f, \"%s\";",DSF_Ave{i,1},DSF_Ave{i,2},DSF_Ave{i,3});
end


fprintf(file_id, "\"%s\", %f, \"%s\"};\n",DSF_Bid{n,1},DSF_Bid{n,2},DSF_Bid{n,3});
fprintf(file_id2, "\"%s\", %f, \"%s\"};\n",DSF_Ask{n,1},DSF_Ask{n,2},DSF_Ask{n,3});
fprintf(file_id3, "\"%s\", %f, \"%s\"};\n",DSF_Ave{n,1},DSF_Ave{n,2},DSF_Ave{n,3});
fclose(file_id);
fclose(file_id2);
fclose(file_id3);
 
endfunction

%Funkcja 8------------------------------------------------------------------------
%funkcja base_curr oblicza walutę BAZOWA
%in: curr_for, curr_dom
%out: base currency
function result =  base_curr(CURR_FOR, CURR_DOM)

order = {"EUR", "GBP", "AUD", "NZD", "USD",  "CAD", "CHF", "NOK", "SEK", "DKK", "CZK", "PLN", "TRY", "MXN", "JPY"};

for i = 1 : columns(order)
if strcmp(order(i), CURR_FOR) == 1
wsk_c = i;
end
if strcmp(order(i), CURR_DOM) == 1
wsk_d = i;
end
end

if wsk_c > wsk_d
tmp = CURR_DOM;
else tmp = CURR_FOR;
end 

result = tmp;

endfunction

%##################Metody interpolacji############################################

%Funkcja 9------------------------------------------------------------------------
%funkcja pomocnicza zwracajaca wektor z wyinterpolowanymi wartosciami
%ulatwia odwolywanie sie do metod interpolacji
%wejscie t-ulamke roku, DS - macierz, 1 kolumna czas, 2 kolumna wartosci DF!
function [dis]=addinterp(t, DS)
global interp_method;

x=DS(:,1);
y=DS(:,2);
n = length(y);

switch (interp_method)
  case "linear on df"
      dis=linear_interp(t,x,y);
  case "linear on rate"
  		kn = (y==1);
  		for i = 1:n
  		if (kn(i) == 1)
  		   y(i) = 0;
  		else  
      	y(i)=df2rate(x(i),y(i));
      endif
      endfor
      temp=linear_interp(t,x,y);
      dis=rate2df(t,temp);
  case "raw"
      dis=raw_interp(t,x,y);
  case "linear on lograte"
      x = x(y<1);
      y = y(y<1); 
      y=log(df2rate(x,y));
      temp=linear_interp(t,x,y);
      dis=rate2df(t,exp(temp));
  case "nat cubic"
     % y=df2rate(x,y); %interpolujemy na stopach proc.
     % temp=ncubs_interp(t,x,y);
     % dis=rate2df(t,temp);
     % interpolacja na discount factors daje lepsze wyniki
     dis = ncubs_interp(t,x,y);
  case "Bspline"
     dis = bspline_interp(t,x,y);   
end

endfunction

%Funkcja 10-----------------------------------------------------------------------
%interpolacja liniowa
%t-wartosci szukane
%(x,y)- wsp. wezla
function [dis]=linear_interp(t,x,y)
t=t(:);%wektor t ma byc pionowy

for i=1:length(t)
      if(t(i)<x(end) && t(i)>x(1)) %interpolacja
	t2=x(x>t(i))(1);
	df2=y(x>t(i))(1);
	t1=x(x<=t(i))(end);
	df1=y(x<=t(i))(end);
	delta=(t(i)-t1)/(t2-t1);
	dis(i)=df2*delta+(1-delta)*df1;	
      elseif (t(i)>=x(end)) %ekstrapolacja, tylko w jedna strone wystarczy bo pierwsza w tablicy start_date
	dis(i)=y(end-1)+((t(i)-x(end-1))/(x(end)-x(end-1)))*(y(end)-y(end-1));
      elseif (t(i)==x(1)) %w przypadku gdy potrzebujemy interpolacji w start_date
	dis(i)=y(1);
      elseif (t(i)<x(1)) %ekstrapolacja w lewa strone
	dis(i)=y(1)-(y(2)-y(1))*(x(1)-t(i))/(x(2)-x(1));
      end
end%for

dis=dis(:);

end

%Funkcja 11-----------------------------------------------------------------------
%interpolacja raw w wersji dla czynnikow dyskontowych  Hagan&West s.96
%t-wartosci szukane
%(x,y)- wsp. wezla
function [dis]=raw_interp(t,x,y)
t=t(:);%wektor t ma byc pionowy

for i=1:length(t)
      if(t(i)<x(end) && t(i)>x(1)) %interpolacja
	t2=x(x>t(i))(1);
	df2=y(x>t(i))(1);
	t1=x(x<=t(i))(end);
	df1=y(x<=t(i))(end);
	delta=(t(i)-t1)/(t2-t1);
	dis(i)=df1^(1-delta)*df2^(delta);	
      elseif (t(i)>=x(end)) %ekstrapolacja, tylko w jedna strone wystarczy bo pierwsza w tablicy start_date
	dis(i)=(y(end))^(t(i)/x(end));
      elseif (t(i)==x(1)) %w przypadku gdy potrzebujemy wartosci w start_date
	dis(i)=y(1);
      end
end%for

dis=dis(:);
end

%Funkcja 12-----------------------------------------------------------------------
%interpolacja natural cubic spline
%algorytm z Kincaid&Cheney Analiza numeryczna s.331/332
%t-wartosci szukane
%(x,y)- wsp. wezla
function [dis]=ncubs_interp(t,x,y)
x=x(:); %zapewnienie wektorow kolumnowych
y=y(:);
t=t(:);

n=length(x);

h=zeros(n,1);
b=zeros(n,1);
u=zeros(n,1);
v=zeros(n,1);

h=diff(x);
b=6*diff(y)./h;


u(1)=2*(h(1)+h(2));
v(1)=b(2)-b(1);

%forward loop
for i=2:n-1;
   u(i)=2*(h(i)+h(i-1))-h(i-1)^2/u(i-1);
   v(i)=b(i)-b(i-1)-h(i-1)*v(i-1)/u(i-1);
end

z(n)=0; 

%backward loop
for i=n-1:-1:2;
   z(i)=(v(i)-h(i)*z(i+1))/u(i);
end
z(1)=0;

k=length(t);

for i=1:k
  if (t(i)<x(end))
    m(i,1)=length(x(x<=t(i))); %indeks okreslajacy poczatek przedzialu
  else
    m(i,1)=n-1; %ostatni przedzial zaczyna sie w t_n-1
  end	%if
end  %for


hm=h(m);
xt1=t-x(m);
xt2=x(m.+ones(k,1))-t;

s1=(z(m)'./(6*hm)).*(xt2).^3;
s2=(z(m.+ones(k,1))'./(6*hm)).*(xt1).^3;
s3=(y(m)./hm.-(z(m)'.*hm)/6).*xt2;
s4=(y(m.+ones(k,1))./hm.-(z(m.+ones(k,1))'.*hm)/6).*xt1;
dis=(s1.+s2.+s3.+s4);

dis=dis(:); %zwracamy wektor kolumnowy

endfunction

%Fu11nkcja 13-----------------------------------------------------------------------
%interpolacja B-spline (generowanie bazy funkcji splinowych)
% na podstawie Lesniewski Interest rates and FX models
function bs = bsp_basis(k,p,t,x) 
% t-argumenty funkcji splinowej, x-wezly na ktorych zbudowano spline, k-stopien spline, 
% p-indeks spline (numer wezla od ktorego spline sie zaczyna) 
bs = zeros(size(t));
if k >= 1
    b1 = bsp_basis(k-1,p,t,x);
    b2 = bsp_basis(k-1,p+1,t,x);

    bs= bs + b1.*((t-x(p))./(x(p+k)-x(p)));
    bs= bs + b2.*((x(p+k+1)-t)./(x(p+k+1)-x(p+1)));

else
    bs(:) = (x(p) <= t).*(t < x(p+1));
end

endfunction


%Funnkcja 14-----------------------------------------------------------------------
%interpolacja B-spline (podstawa krzywej splinowej)% algorytm na podstawie James&Webber Interest Rate Modelling str. 434-438 z poprawkami we wzorze (15.26)
%t-wartosci szukane
%(x,y)- wsp. wezla 

function [dis]=bspline_interp(t,x,y,deg)
%deg = stopien
if nargin<4
 deg=3;
end

z=[-5,-3,-1,0,1,3,5,7,10,15,20,30,40,50]; %siatka punktow generacji spline. Tutaj uwazac bo x w srodku ma byc wekorem jednowierszowym

knt = length(z)-4;

B = zeros(length(x),knt);
for j = 1 : knt
    B(:,j) = bsp_basis(deg,j,x,z);
end

%c=(B'*B)\(B'*y);
c=pinv(B)*y;
%c=B\y;

W = zeros(knt, 1);
for j = 1 : knt
    W(j) = bsp_basis(deg,j,0,z);
end

tm = pinv(B'*B) * W;
cs = c + (1 - W' *c)/(W'* tm)* tm;

%cs = c + ((B'*B)\ W)/( W' * ((B'*B)\ W))*(1 - W' *c);

D=zeros(length(t),knt);
for j = 1 : knt
    D(:,j) = bsp_basis(deg,j,t,z);
end

dis = D * cs;
dis=dis(:); %zwracamy wektor kolumnowy

endfunction




%############Funkcje liczące DF i FR##############################################

%Funkcja 15-----------------------------------------------------------------------
%czynnik dyskontowy na chwilę bieżącą
%t - czas ulamek roku
%DS - tablica czynnikow DF
function [dis]=DF_pom(t,DS)
global start_date;

n=length(DS);
for i=1:n;
  DS1(i, 1)=year_frac(start_date, DS{i, 1}, DS{i,3});
  DS1(i, 2)=DS{i, 2};
end


dis=addinterp(t,DS1);

endfunction


%Funkcja 16-----------------------------------------------------------------------
%czynnik dyskontowy forward
function [y]=DF(date0, date,DS)
global start_date;

T0=year_frac(start_date, date0, DS{DS_conv(date0,DS),3});
n=length(date);
y=zeros(n, 1);
	for i=1:n;
		idx=DS_conv(date{i},DS);
		a = DF_pom(year_frac(start_date, date{i}, DS{idx,3}),DS);
		if (T0 == 0)
			b = 1;
		else
			b = DF_pom(T0,DS);
		endif
		y(i)= a/b;
	endfor;
y=y';
endfunction

%Funkcja 17-----------------------------------------------------------------------
%stopa forward
function [x]=FR(date0, dates1, dates2,DS)
df1=DF(date0, dates1,DS);
df2=DF(date0, dates2,DS);

n=length(dates1);
	for i=1:n
		idx=DS_conv(dates1{i},DS);
		idx2=DS_conv(dates2{i},DS);
		time(i)=year_frac(date0, dates2{i}, DS(idx2,3))-year_frac(date0, dates1{i}, DS(idx,3));
	endfor
x=(1./time).*(df1./df2-1);
endfunction


%Funkcja 18-----------------------------------------------------------------------
%funkcja pomocnicza zwracająca konwencje
function [idx]=DS_conv(date, DS)
  for i=1:length(DS)
    idx=i;
    if (day_diff(date, DS{i,1},"ACT")>=0)
      break;
    end
  end

endfunction


%Funkcja 19 ---------------------------------------------------------------------
% konwertuje tablice DOM i FOR do BAS i QUO

function convert_tables()
 
global CURR_DOM;
global CURR_FOR;
global BAS_CURR;
global QUO_CURR;
 
global DSD_Bid;
global DSD_Ask;
global DSD_Ave;
global DSF_Ask;
global DSF_Bid;
global DSF_Ave;

global DSQ_Ask;
global DSQ_Bid;
global DSQ_Ave;
global DSB_Ask;
global DSB_Bid;
global DSB_Ave; 
 
source DSD_Bid.m;
source DSD_Ask.m;
source DSD_Ave.m;
source DSF_Ask.m;
source DSF_Bid.m;
source DSF_Ave.m;

if (BAS_CURR==CURR_DOM) %jesli bazowa DOM
     DSQ_Bid = DSF_Bid;  
     DSQ_Ask = DSF_Ask;  
     DSQ_Ave = DSF_Ave;   
     DSB_Bid = DSD_Bid;  
     DSB_Ask = DSD_Ask;  
     DSB_Ave = DSD_Ave;   
elseif(BAS_CURR==CURR_FOR) %jeśli bazowa jest FOR
     DSQ_Bid = DSD_Bid;  
     DSQ_Ask = DSD_Ask;  
     DSQ_Ave = DSD_Ave;  
     DSB_Bid = DSF_Bid;  
     DSB_Ask = DSF_Ask;  
     DSB_Ave = DSF_Ave;   
end   

n=length(DSQ_Bid);

file_id = fopen('DSQ_Bid.m', 'w');
file_id2 = fopen('DSQ_Ask.m', 'w');
file_id3 = fopen('DSQ_Ave.m', 'w');

fprintf(file_id, "DSQ_Bid={");
fprintf(file_id2, "DSQ_Ask={");
fprintf(file_id3, "DSQ_Ave={");

for i=1:n-1
    fprintf(file_id, "\"%s\", %f, \"%s\";",DSQ_Bid{i,1},DSQ_Bid{i,2},DSQ_Bid{i,3});
    fprintf(file_id2, "\"%s\", %f, \"%s\";",DSQ_Ask{i,1},DSQ_Ask{i,2},DSQ_Ask{i,3});
    fprintf(file_id3, "\"%s\", %f, \"%s\";",DSQ_Ave{i,1},DSQ_Ave{i,2},DSQ_Ave{i,3});
end


fprintf(file_id, "\"%s\", %f, \"%s\"};\n",DSQ_Bid{n,1},DSQ_Bid{n,2},DSQ_Bid{n,3});
fprintf(file_id2, "\"%s\", %f, \"%s\"};\n",DSQ_Ask{n,1},DSQ_Ask{n,2},DSQ_Ask{n,3});
fprintf(file_id3, "\"%s\", %f, \"%s\"};\n",DSQ_Ave{n,1},DSQ_Ave{n,2},DSQ_Ave{n,3});
fclose(file_id);
fclose(file_id2);
fclose(file_id3);

n=length(DSB_Bid);

file_id = fopen('DSB_Bid.m', 'w');
file_id2 = fopen('DSB_Ask.m', 'w');
file_id3 = fopen('DSB_Ave.m', 'w');

fprintf(file_id, "DSB_Bid={");
fprintf(file_id2, "DSB_Ask={");
fprintf(file_id3, "DSB_Ave={");

for i=1:n-1
    fprintf(file_id, "\"%s\", %f, \"%s\";",DSB_Bid{i,1},DSB_Bid{i,2},DSB_Bid{i,3});
    fprintf(file_id2, "\"%s\", %f, \"%s\";",DSB_Ask{i,1},DSB_Ask{i,2},DSB_Ask{i,3});
    fprintf(file_id3, "\"%s\", %f, \"%s\";",DSB_Ave{i,1},DSB_Ave{i,2},DSB_Ave{i,3});
end


fprintf(file_id, "\"%s\", %f, \"%s\"};\n",DSB_Bid{n,1},DSB_Bid{n,2},DSB_Bid{n,3});
fprintf(file_id2, "\"%s\", %f, \"%s\"};\n",DSB_Ask{n,1},DSB_Ask{n,2},DSB_Ask{n,3});
fprintf(file_id3, "\"%s\", %f, \"%s\"};\n",DSB_Ave{n,1},DSB_Ave{n,2},DSB_Ave{n,3});
fclose(file_id);
fclose(file_id2);
fclose(file_id3);



endfunction



 
 

%-----------------------------------------------------

% funkcja 20 (wyznaczanie trzeciej srody miesiaca)
% in: data1 - data pierwszego dnia miesiaca one_tenor, one_tenor (miesiac w ktorym rozpoczyna sie dany kontrakt futures)
% out: data trzeciej srody danego miesiaca, nazwa miesiaca  

function [data,month] = wedn(date1,one_tenor)
	i=weekday(date1);
	if (i<5)
		data=num2str(19-i);
	else
		data=num2str(26-i);
	end
	
	if (strncmp(one_tenor,"Dec",3)==1)
		month="Dec";
	elseif (strncmp(one_tenor,"Mar",3)==1)
		month="Mar";
	elseif (strncmp(one_tenor,"Jun",3)==1)
		month="Jun";
	elseif (strncmp(one_tenor,"Sep",3)==1)
		month="Sep";
	end
		
end 