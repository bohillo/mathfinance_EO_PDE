1;
% FINANSE OBLICZENIOWE
% DUŻY PROJEKT 2010 i 2011 - FUNKCJE KALENDARZOWE
% A deadline is negative inspiration. Still, it's better than no inspiration at all.  - Rita Mae Brown 
%
% WERSJA 1.02 (final, known bugs fixed)
% 1 day = 60*60*24 sec = 86400 sec
% 1 hour = 60*60 sec = 3600 sec
% WERSJA 2.0 - DUZY PROJEKT 2011

% funkcja 1
% in: ustalona data w formacie dd-mmm-yyyy np. 07-Sep-2000
% out: dzien tygodnia
% UWAGA: juz zaimplementowana w Octave: weekday("dd-mmm-yyyy") (od: 1=niedziela do: 7=sobota)

% funkcja 2
% in: data1 (poczatkowa), data2 (koncowa), data1<=data2, formatu jak wyzej, konwencja dni
% out: liczba dni pomiedzy datami w sensie dane konwencji
%
% identyfikatory konwencji dni: ACT, 30, 30E, 30E+

function days = day_diff(date1,date2,basis)
	tm1 = strptime(date1, "%d-%b-%Y");  
	tm2 = strptime(date2, "%d-%b-%Y");
	
	if strcmp(basis,"ACT") 
	days = round((mktime(tm2)-mktime(tm1))/86400);
	
	elseif strcmp(basis,"30") 
	if(tm1.mon==1 && tm1.mday==eomday(tm1.year+1900,2)) tm1.mday=30; if(tm2.mon==1 && tm2.mday==eomday(tm2.year+1900,2)) tm2.mday=30; end end
	if(tm1.mday==31) tm1.mday=30; end
	if(tm2.mday==31 && (tm1.mday==30 ||tm1.mday==30)) tm1.mday=30; end 
	days = (tm2.year-tm1.year)*360+(tm2.mon-tm1.mon)*30+tm2.mday-tm1.mday;
	
	elseif strcmp(basis,"30E") 
	days = (tm2.year-tm1.year)*360+(tm2.mon-tm1.mon)*30+min(tm2.mday,30)-min(tm1.mday,30);
	
	elseif strcmp(basis,"30E+")
	if(tm2.mday==31) tm2.mday=1; tm2.mon=mod(tm2.mon+1,11); end 
	days = (tm2.year-tm1.year)*360+(tm2.mon-tm1.mon)*30+tm2.mday-min(tm1.mday,30);
	
	% ewentualne dalsze konwencje
	end
end

% funkcja 3
% in: data1, formatu jak wyżej, liczba dni n
% out: data uzyskana przez operacje: data1 + n dni

function date2 = add_days(date1,n)
	tm1 = strptime(date1, "%d-%b-%Y");
	date2 = strftime("%d-%b-%Y",localtime(mktime(tm1)+n*86400+7200));
end

% funkcja 4
% in: data1, wektor liczby miesiecy w okresach (1Y = 12M), parametr End of Month Adjustment +/-1 ("w przod"/"w tyl")
% out: lista dat wyznaczonych z danych wejsciowych

function dat_list = add_months(date1,mwect,EMA)
	months={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
	n = length(mwect);
	w = cumsum(mwect);
	tm1 = strptime(date1, "%d-%b-%Y");
	y = tm1.year; m = tm1.mon; d = tm1.mday;
	dat_list = cell(n+1,1);
	dat_list{1} = date1;
	for i=2:n+1
		tm1.mon = mod(m+w(i-1),12);
		tm1.year = y + floor((m+w(i-1))/12);
		t=eomday(tm1.year+1900,tm1.mon+1);
		if (t<d)
			if (EMA==-1) tm1.mday=t;
			elseif (EMA==1) tm1.mday=1; tm1.year=tm1.year+(tm1.mon==11); tm1.mon=mod(tm1.mon+1,12); 
			end
		else tm1.mday = d;	
		end
		day=["0" int2str(tm1.mday)];
		dat_list{i} = [ day([length(day)-1,length(day)]) "-" months{tm1.mon+1} "-" int2str(tm1.year+1900)];
	end
end

% funkcja 5 (zmieniona)
% in: rok, identyfikator FC
% out: lista bank holidays wyznaczonych przez Wielkanoc na danej gieldzie

function easter_hol = easter(year,FC)
	FC = lower(FC);
	
	%implementacja algorytmu Meeusa/Jonesa/Butchera
	a = mod(year,19);
	b = floor(year/100);
	c = mod(year,100);
	h = mod(19*a+b-b/4-floor((b-floor((b+8)/25)+1)/3)+15,30);
	l = mod(32+2*mod(b,4)+2*floor(c/4)-h-mod(c,4),7);
	m = floor((a+11*h+22*l)/451);
	d = mod(h+l-7*m+114,31)+1; %dzien 
	m = floor((h+1-7*m+114)/31); %miesiac
	% najwczesniej 22 marca, najpozniej 25 kwietnia

	tm1 = strptime("01-Jan-2000", "%d-%b-%Y");
	tm1.mday = d; tm1.mon = m-1; tm1.year = year-1900;
	
	if strcmp("warsaw",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+60*86400+7200))};
	% Wlk. Piatek, Pon. Wielk. i Boze Cialo
	
	elseif strcmp("copenhagen",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-3*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)),strftime("%d-%b-%Y",localtime(mktime(tm1)+26*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+39*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+40*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+50*86400+7200))};
	% Maundy Thursday, Good Friday, Easter Monday, General Prayer Day, Ascension Thursday, Friday after Ascension Thursday(Bank Holiday), Whit Monday 
	
	elseif strcmp("london",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200))};
	% Good Friday, Easter Monday
	
	elseif strcmp("zurich",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+39*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+50*86400+7200))};
	% Good Friday, Easter Monday, Ascension Day, Whit Monday
	
	elseif strcmp("frankfurt",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+39*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+50*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+60*86400+7200))};
	% Good Friday, Easter Monday, Ascension Thursday, Whit Monday, Corpus Christi
	
	elseif strcmp("oslo",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-3*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+39*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+50*86400+7200))};
	% Holy Thursday, Good Friday, Easter Monday, Ascension Thursday, Whit Monday 
	
	elseif strcmp("stockholm",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+39*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+50*86400+7200))};
	% Good Friday, Easter Monday, Ascension Thursday, Whit Monday 
	
	elseif strcmp("target",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200))};
	% Good Friday, Easter Monday
	
	elseif strcmp("sydney",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200))};
	% Good Friday, Easter Monday
	
	elseif strcmp("toronto",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200))};
	% Good Friday, Easter Monday
	
	elseif (strcmp("newyork",FC)| strcmp("nymex",FC))
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200))};
	% Good Friday
	
	elseif strcmp("wellington",FC)
	easter_hol = {strftime("%d-%b-%Y",localtime(mktime(tm1)-2*86400+7200)), strftime("%d-%b-%Y",localtime(mktime(tm1)+86400+7200))};
	% Good Friday, Easter Monday
	end
	
end

% funkcja 6 (zmieniona)
% in: data, identyfikator FC
% out: 0 lub 1 (dany dzien nie jest/jest dniem roboczym)

function is = is_business_day(date,FC)
	FC = lower(FC);
	tm = strptime(date, "%d-%b-%Y");
	%obliczenia rownonocy wiosennej i jesiennej - potrzebne dla dni wolnych w Tokio:
	%zrodlo: QuantLib
	  exact_vernal_equinox_time = 20.69115;
	  exact_autumnal_equinox_time = 23.09;
	  moving_amount = (tm.year-100)*0.242194;
          number_of_leap_years = floor((tm.year-100)/4+(tm.year-100)/100-(tm.year-100)/400);
          ve = floor(exact_vernal_equinox_time + moving_amount - number_of_leap_years);
          ae = floor(exact_autumnal_equinox_time + moving_amount - number_of_leap_years);
	  if (tm.year+1900)==1997 ae=22; endif % wyjatek wynikajacy z niedokladnosci algorytmu
	
	if strcmp("warsaw",FC)
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% Nowy Rok
	|| (tm.mday==6 && tm.mon==0 && (tm.year+1900)>=2011)                    % Trzech Kroli - od 2011
	|| (tm.mday==1 && tm.mon==4)						% sw. Panstwowe
	|| (tm.mday==3 && tm.mon==4)						% sw. Konst. 3 Maja
	|| (tm.mday==15 && tm.mon==7)						% Wnieb. NPM
	|| (tm.mday==1 && tm.mon==10)						% Wsz. swietych
	|| (tm.mday==11 && tm.mon==10)						% sw. Niepodleglosci
	|| (tm.mday==24 && tm.mon==11)						% Wigilia
	|| (tm.mday==25 && tm.mon==11)						% Boze Narodzenie I
	|| (tm.mday==26 && tm.mon==11)						% Boze Narodzenie II
	|| ismember(date,easter(tm.year+1900,"warsaw")));	% Wlk. Piatek, Pon. Wielk. i Boze Cialo

   elseif strcmp("copenhagen",FC)
	is = ~(ismember(weekday(date),[1,7])				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==5 && tm.mon==5)						% Constitution Day
	|| (tm.mday==24 && tm.mon==11)						% Christmas Eve
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11)						% Christmas
	|| (tm.mday==31 && tm.mon==11)						% New Year's Eve
	|| ismember(date,easter(tm.year+1900,"copenhagen")));	% Maundy Thursday, Good Friday, Easter Monday, General Prayer Day, Ascension Day,  Friday after Ascension Thursday(Bank Holiday), Whit Monday

	elseif strcmp("london",FC)
	is = ~(ismember(weekday(date),[1,7])						% niedziela, sobota
	|| ((tm.mday==1 || ((tm.mday== 2 || tm.mday==3) && tm.wday==1)) && tm.mon==0)		% New Year's Day 
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==4)						% Early May Bank Holiday
	|| (tm.mday>=25 && tm.wday==1 && tm.mon==4 && (tm.year+1900)~= 2002)						% Spring Bank Holiday
	|| (tm.mday>=25 && tm.wday==1 && tm.mon==7 )						% Summer Bank Holiday
	|| ((tm.mday==25 || (tm.mday==27 && (tm.wday==1 || tm.wday==2))) && tm.mon==11)		% Christmas Day
	|| ((tm.mday==26 || (tm.mday==28 && (tm.wday==1 || tm.wday==2))) && tm.mon==11)	 	% Boxing Day
	|| ((tm.mday==3 || tm.mday==4) && tm.mon==5 && (tm.year+1900)==2002)  % Golden Jubilee Bank Holiday & special Spring Bank Holiday
	|| (tm.mday==31 && tm.mon==11 && (tm.year+1900)==1999) % December 31st, 1999 only
	|| ismember(date,easter(tm.year+1900,"london")));				% Good Friday, Easter Monday
	
	elseif strcmp("zurich",FC)
	is = ~(ismember(weekday(date),[1,7])				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==2 && tm.mon==0) 						% Berchtoldstag
	|| (tm.mday==1 && tm.mon==4)						% Labour Day
	|| (tm.mday==1 && tm.mon==7)						% National Day
	|| (tm.mday==24 && tm.mon==11)						% Christmas Eve
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11)						% St. Stephen's Day
	|| (tm.mday==31 && tm.mon==11)						% New Year's Eve
	|| ismember(date,easter(tm.year+1900,"zurich")));	% Good Friday, Easter Monday, Ascension Day, Whit Monday
		
	elseif strcmp("frankfurt",FC)
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==1 && tm.mon==4)						% Labour Day
	|| (tm.mday==3 && tm.mon==9)						% National Day
	|| (tm.mday==24 && tm.mon==11)						% Christmas Eve
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11)						% Boxing Day
	|| (tm.mday==31 && tm.mon==11)						% New Year's Eve
	|| ismember(date,easter(tm.year+1900,"frankfurt")));% Good Friday, Easter Monday, Ascension Thursday, Whit Monday, Corpus Christi
	
	elseif (strcmp("newyork",FC)| strcmp("nymex",FC))
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==2 && tm.mon==0 && tm.wday==1)				% New Year's Day moved to Monday if on Sunday
	|| (tm.mday==31 && tm.mon==11 && tm.wday==5)				% New Year's Day moved to Friday if on Saturday
	|| (tm.mday>=15 && tm.mday<=21 && tm.mon==0 && tm.wday==1 && (tm.year+1900)>=1998 )		% Martin Luther King's birthday
	|| (tm.mday>=15 && tm.mday<=21 && tm.mon==1 && tm.wday==1)		% Washington's birthday
	|| (tm.mday>=25 && tm.mon==4 && tm.wday==1)				% Memorial Day
	|| (tm.mday==4 && tm.mon==5)						% Independence Day
	|| (tm.mday==5 && tm.mon==5 && tm.wday==1)				% Independence Day moved to Monday if on Sunday
	|| (tm.mday==3 && tm.mon==5 && tm.wday==5)				% Independence Day moved to Friday if on Saturday
	|| (tm.mday<=7 && tm.mon==8 && tm.wday==1)				% Labour Day
	|| (tm.mday>=22 && tm.mday<=28 && tm.mon==10 && tm.wday==4)		% Thanksgiving Day
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11 && tm.wday==1)				% Christmas moved to Monday if on Sunday
	|| (tm.mday==24 && tm.mon==11 && tm.wday==5)				% Christmas moved to Friday if on Saturday
	|| (tm.mday==11 && tm.mon==5 && (tm.year+1900)==2004 ) % Reagan's funeral
	|| (tm.mday==11 && tm.mon==8 && (tm.year+1900)==2001 ) % September 11, 2001
	|| (tm.mday==14 && tm.mon==6 && (tm.year+1900)==1977 ) % 1977 Blackout
	|| (tm.mday==25 && tm.mon==0 && (tm.year+1900)==1973 ) % Lyndon B. Johnson funeral
	|| (tm.mday==28 && tm.mon==11 && (tm.year+1900)==1972 ) % Harry S. Truman funeral
	|| ismember(date,easter(tm.year+1900,"newyork")));% Good Friday

   elseif strcmp("oslo",FC)
	is = ~(ismember(weekday(date),[1,7])				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==1 && tm.mon==4)						% Labour Day
	|| (tm.mday==17 && tm.mon==4)						% National Day
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11)						% Boxing Day
	|| ismember(date,easter(tm.year+1900,"oslo")));	% Holy Thursday, Good Friday, Easter Monday, Ascension Day, Whit Monday

   elseif strcmp("stockholm",FC)
	is = ~(ismember(weekday(date),[1,7])				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==6 && tm.mon==0)                 % Epiphany 
	|| (tm.mday==1 && tm.mon==4)						% May Day
	|| (tm.mday==6 && tm.mon==5 && (tm.year+1900)>=2011) % National Day since 2011
	|| (tm.mday>=18 && tm.mday<=24 && tm.mon==5 && tm.wday==5)					% Midsummer Eve
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11)						% Boxing Day
	|| (tm.mday==31 && tm.mon==11)						% New Year's Eve
	|| ismember(date,easter(tm.year+1900,"stockholm")));	% Good Friday, Easter Monday, Ascension Day, Whit Monday


	elseif strcmp("sydney",FC)
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0)						% New Year's Day
	|| ((tm.mday==2 || tm.mday==3) && (tm.mon==0 && tm.wday==1))		% New Year's Day moved to Monday
	|| (tm.mday==26 && tm.mon==0)						% Australia Day
	|| ((tm.mday==27 || tm.mday==28) && (tm.mon==0 && tm.wday==1))		% Australia Day moved to Monday
	|| (tm.mday==25 && tm.mon==3)						% ANZAC Day
	|| (tm.mday==26 && tm.mon==3 && tm.wday==1)				% ANZAC Day moved to Monday
	|| (tm.mday>7 && tm.mday<=14 && tm.wday==1 && tm.mon==5)		% Queen's Birthday
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==7)				% Bank Holiday
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==9)				% Labour Day
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==27 && tm.mon==11 && (tm.wday==1 || tm.wday==2))		% Christmas moved to Monday or Tuesday
	|| (tm.mday==26 && tm.mon==11)						% Boxing Day
	|| (tm.mday==28 && tm.mon==11 && (tm.wday==1 || tm.wday==2))		% Boxing Day moved to Monday or Tuesday
	|| ismember(date,easter(tm.year+1900,"sydney")));	% Good Friday, Easter Monday 
	
	elseif strcmp("target",FC)
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0)						% New Year's Day
	|| (tm.mday==1 && tm.mon==4 && (tm.year+1900)>=2000)			% Labour Day - od 2000
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==26 && tm.mon==11 && (tm.year+1900)>=2000)			% Day of Goodwill - od 2000
	|| (tm.mday==31 && tm.mon==11 && (tm.year+1900==1998 || tm.year+1900==1999 || tm.year+1900==2000)) %31 grudnia tylko w tych latach
	|| ismember(date,easter(tm.year+1900,"target")));	% Good Friday, Easter Monday 
          
	elseif strcmp("tokyo",FC) 
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0)						% New Year's Day
	|| (tm.mday==2 && tm.mon==0)						% Bank Holiday
	|| (tm.mday==3 && tm.mon==0)						% Bank Holiday
	|| (tm.mday>=8 && tm.mday<=14 &&tm.wday==1 && (tm.year+1900)>=2000 && tm.mon==0)	% Coming of Age Day
	|| ((tm.mday==15 || (tm.mday==16 && tm.wday==1)) && tm.mon==0 && (tm.year+1900)<2000)	% Coming of Age Day before 2000
	|| ((tm.mday==11 || (tm.mday==12 && tm.wday==1)) && tm.mon==1)		% National Foundation Day
	|| ((tm.mday==ve || (tm.mday==ve+1 && tm.wday==1)) && tm.mon==2)	% Vernal Equinox
	|| ((tm.mday==29 || (tm.mday==30 && tm.wday==1)) && tm.mon==3)		% Showa Day
	|| (tm.mday==3 && tm.mon==4)						% Constitution Memorial Day
	|| (tm.mday==4 && tm.mon==4)						% Greenery Day
	|| ((tm.mday==5 || (tm.mday==6 && tm.wday==1)) && tm.mon==4)		% Children's Day
	|| (tm.mday>=15 && tm.mday<=21 && tm.wday==1 && tm.mon==6 && (tm.year+1900)>=2003) % Marine Day
	|| ((tm.mday==20 || (tm.mday==21 && tm.wday==1)) && tm.mon==6 && (tm.year+1900)>=1996 && (tm.year+1900)<2003)	%Marine Day between 1996 and 2003
	|| (tm.mday>=15 && tm.mday<=21 && tm.wday==1 && tm.mon==8 && (tm.year+1900)>=2003) % Respect for the Aged Day
	|| ((tm.mday==15 || (tm.mday==16 && tm.wday==1)) && tm.mon==8 && (tm.year+1900)<2003) % Respect for the Aged Day before 2003
	|| (tm.mday==ae-1 && tm.wday==2 && tm.mday>=16 && tm.mday<=22 && tm.mon==8 && (tm.year+1900)>=2003) % Exactly one day between Respect for the Aged Day and Autumnal Equinox
	|| ((tm.mday==ae || (tm.mday==ae+1 && tm.wday==1)) && tm.mon==8 )	% Autumnal Equinox
	|| (tm.mday>=8 && tm.mday <=14 && tm.wday==1 && tm.mon==9 && (tm.year+1900)>=2000) % Health and Sports Day
	|| ((tm.mday==10 || (tm.mday==11 && tm.wday==1)) && tm.mon==9 && (tm.year+1900)<2000) % Health and Sports Day before 2000
	|| ((tm.mday==3 || (tm.mday==4 && tm.wday==1)) && tm.mon==10)		% Culture Day
	|| ((tm.mday==23 || (tm.mday==24 && tm.wday==1)) && tm.mon==10)		% Labor Thanksgiving Day
	|| ((tm.mday==23 || (tm.mday==24 && tm.wday==1)) && tm.mon==11 && (tm.year+1900)>=1989) % Emperor's Birthday
	|| (tm.mday==31 && tm.mon==11)						% Bank Holiday
	|| (tm.mday==10 && tm.mon==3 && (tm.year+1900)==1959)			% Marriage of Prince Akihito
	|| (tm.mday==24 && tm.mon==1 && (tm.year+1900)==1989)			% Rites of Imperial Funeral
	|| (tm.mday==12 && tm.mon==10 && (tm.year+1900)==1990)			% Enthronement Ceremony
	|| (tm.mday==9 && tm.mon==5 && (tm.year+1900)==1993));			% Marriage of Prince Naruhito	           
	
	elseif strcmp("toronto",FC) 
	is = ~(ismember(weekday(date),[1,7]) 				% niedziela, sobota
	|| (tm.mday==1 && tm.mon==0) 						% New Year's Day
	|| (tm.mday==2 && tm.mon==0 && tm.wday==1)				% New Year's Day moved to Monday
	|| (tm.mday>17 && tm.mday<=24 && tm.wday==1 && tm.mon==4)		% Victoria Day
	|| (tm.mday==1 && tm.mon==6)						% Canada Day
	|| ((tm.mday==2 || tm.mday==3) && (tm.mon==6 && tm.wday==1))		% Canada Day moved to Monday
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==7)				% Provincial Holiday
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==8)				% Labor Day
	|| (tm.mday>7 && tm.mday<=14 && tm.wday==1 && tm.mon==9)		% Thanksgiving Day
	|| (tm.mday==11 && tm.mon==10)						% November 11th
	|| (tm.mday==25 && tm.mon==11)						% Christmas
	|| (tm.mday==27 && tm.mon==11 && (tm.wday==1 || tm.wday==2))		% Christmas moved to Monday or Tuesday
	|| (tm.mday==26 && tm.mon==11)						% Boxing Day
	|| (tm.mday==28 && tm.mon==11 && (tm.wday==1 || tm.wday==2)) 		% Boxing Day moved to Monday or Tuesday
	|| ismember(date,easter(tm.year+1900,"toronto")));				% Good Friday, Easter Monday
	
	elseif strcmp("wellington",FC)
	is = ~(ismember(weekday(date),[1,7])						% niedziela, sobota
	|| ((tm.mday==1 || (tm.mday==3 && (tm.wday==1 || tm.wday==2))) && tm.mon==0)		% New Year's Day 
	|| ((tm.mday==2 || (tm.mday==4 && (tm.wday==1 || tm.wday==2))) && tm.mon==0)		% Day after New Year's Day 
	|| (tm.mday>=19 && tm.mday<=25 && tm.wday==1 && tm.mon==0)						% Anniversary Day
	|| (tm.mday==6 && tm.mon==1)						% Waitangi Day
	|| (tm.mday==25 && tm.mon==3)						% ANZAC Day
	|| (tm.mday<=7 && tm.wday==1 && tm.mon==5)	% Queen's Birthday
	|| (tm.mday>=22 && tm.mday<=28 && tm.wday==1 && tm.mon==9)  % Labour Day
	|| ((tm.mday==25 || (tm.mday==27 && (tm.wday==1 || tm.wday==2))) && tm.mon==11)		% Christmas Day
	|| ((tm.mday==26 || (tm.mday==28 && (tm.wday==1 || tm.wday==2))) && tm.mon==11)	 	% Boxing Day
	|| ismember(date,easter(tm.year+1900,"wellington")));	% Good Friday, Easter Monday
	
	
	% ewentualne dalsze giełdy
	end
end

% funkcja 6.2
% in: data, cell array identyfikatorow FC
% out: 0 lub 1 (dane dni nie sa/sa robocze (wszystkie, co najmniej jeden, wybrane kilka z nich - w zaleznosci od postaci FCL))
function is = is_business_day2(date,FCL)
  if iscell(FCL{length(FCL)})
    is=1;
    n=length(FCL);
    if FCL{n}{1}==0
      for i=1:(n-1) 
	is= is || is_business_day(date,FCL{i});
	endfor;
    else
      for i=1:length(FCL{n}) 
	is= is && is_business_day(date,FCL{FCL{n}{i}});
	endfor;  
    endif;
  else
    is=1;
    for i=1:length(FCL)
      is= is && is_business_day(date,FCL{i});
      endfor;
  endif;
end


% funkcja 7 (zmieniona)
% in: data, cell array identyfikatorow FC, parametr +/-1 ("w przod"/"w tyl")
% out: nablizszy nastepny/poprzedni dzien roboczy 

function date2 = day_shift(date1,FCL,x)
	date2 = date1;
	do date2 = add_days(date2,x); until is_business_day2(date2,FCL);
end

% funkcja 8 (modyfikacja 7) (zmieniona)
% in: data1, cell array identyfikatorow FC, liczba dni (roboczych)
% out: data2 powstala przez przesuniecie daty1 o n dni roboczych do przodu / do tylu

function date2 = day_shift2(date1,FCL,n)
	date2 = date1;
   if (n>0)
   	k = 0;
	   do date2 = add_days(date2,1); k = k + is_business_day2(date2,FCL); until (k==n);
	end 
	if (n<0)
   	k = 0;
	   do date2 = add_days(date2,-1); k = k - is_business_day2(date2,FCL); until (k==n);
	end 
end

% funkcja 9 (zmieniona)
% in: data, cell array identyfikatorow FC ,parametr +/-1
% out: nablizszy nastepny/poprzedni dzien roboczy w sensie konwencji Modified Following/Previous Business Day

function date2 = mod_day_shift(date1,FCL,x)
	date2 = day_shift(date1,FCL,x);
	tm1 = strptime(date1, "%d-%b-%Y");
	tm2 = strptime(date2, "%d-%b-%Y");
	if (tm1.mon~=tm2.mon) date2 = day_shift(date1,FCL,-x); end
end
	
% funkcja 10 (zmieniona)
% in: data1, cell array identyfikatorow FC, konwencja Business Day Adjustment
% out: dzien platnosci odpowiadajacy data1 (wyznaczony przez konwencje) tzn. odpowiednia data2
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual

function date2 = date_rolling(date1,FCL,BDA)
	if is_business_day2(date1,FCL) 
		date2 = date1;
	else 
		if  strcmp(BDA,"actu")
			date2 = date1;
			
		elseif strcmp(BDA,"sfbd")
			date2 = day_shift(date1,FCL,1);
			
		elseif strcmp(BDA,"mfbd") 
			date2 = mod_day_shift(date1,FCL,1);
			
		elseif strcmp(BDA,"spbd") 
			date2 = day_shift(date1,FCL,-1);
			
		elseif strcmp(BDA,"mpbd")
			date2 = mod_day_shift(date1,FCL,-1);
		
		elseif strcmp(BDA,"seom") 
			tm1 = strptime(date1, "%d-%b-%Y");
			tm1.mday = eomday(tm1.year+1900,tm1.mon+1);
			date2 = strftime("%d-%b-%Y",localtime(mktime(tm1)));
			date2 = date_rolling(date2,FCL,"spbd");
		
		% ewentualne dalsze konwencje
		end	
	end
end

% funkcja 11
% in: data1, data2, formatu jak wyzej
% out: ulamek roku (w konwencji ACT/ACT) wyznaczony przez te daty 

function fraction = act_act_frac(date1,date2)
	tm1 = strptime(date1, "%d-%b-%Y");
	tm2 = strptime(date2, "%d-%b-%Y");
	R1=tm1.year+1900; R2=tm2.year+1900;
	
	if (R1==R2) 
	fraction = day_diff(date1,date2,"ACT")/(365+is_leap_year(R1));
	else 
	fraction = day_diff(date1,["31-Dec-" num2str(R1)],"ACT")/(365+is_leap_year(R1)) + R2-R1-1 + day_diff(["31-Dec-" num2str(R2-1)],date2,"ACT")/(365+is_leap_year(R2));
	end
end

% funkcja 11.1
% in: data1, data2, formatu jak wyzej
% out: ulamek roku (w konwencji ACT/ACT AFB) wyznaczony przez te daty 
function fraction = act_act_afb_frac(date1,date2)
	tm1 = strptime(date1, "%d-%b-%Y");
	tm2 = strptime(date2, "%d-%b-%Y");
	R1=tm1.year+1900; R2=tm2.year+1900;
	
	if (R1==R2 || (R1+1==R2 && (tm1.mon>tm2.mon || (tm1.mon==tm2.mon && tm1.mday>=tm2.mday))))
	  % odcinek czasu max 1 rok i kolejno: zawiera/nie zawiera 29.02 w przedziale (date1,date2]
	  if ((is_leap_year(R1) && (tm1.mon==0 || (tm1.mon==1 && tm1.mday<29)) && (tm2.mon>=2 || R2>R1)) 
	  || (is_leap_year(R2) && (tm2.mon>1 || (tm2.mon==1 && tm2.mday==29)) && ((tm1.mon==0 || (tm1.mon==1 && tm1.mday<29)) || R2>R1)))
	  fraction = day_diff(date1,date2,"ACT")/366;
	  else
	  fraction = day_diff(date1,date2,"ACT")/365;
	  end
	else
	  % odcinek czasu dluzszy niz 1 rok i j/w
	  fraction = 0;
	  do
	    fraction = fraction + 1;
	    date2=strrep(date2,int2str(R2),int2str(R2-1));
	    R2 = R2 - 1;
	  until (R1==R2 || (R1+1==R2 && (tm1.mon>tm2.mon || (tm1.mon==tm2.mon && tm1.mday>=tm2.mday))))
	  if (is_leap_year(R2) && tm2.mon==1 && tm2.mday==28) tm2.mday==29; date2=strrep(date2,"28-", "29-"); end;
	  if (!is_leap_year(R2) && tm2.mon==1 && tm2.mday==29) tm2.mday==28; date2=strrep(date2,"29-", "28-"); end;
	  
	  if ((is_leap_year(R1) && (tm1.mon==0 || (tm1.mon==1 && tm1.mday<29)) && (tm2.mon>=2 || R2>R1)) 
	  || (is_leap_year(R2) && (tm2.mon>1 || (tm2.mon==1 && tm2.mday==29)) && ((tm1.mon==0 || (tm1.mon==1 && tm1.mday<29)) || R2>R1)))
	  fraction = fraction + day_diff(date1,date2,"ACT")/366;
	  else
	  fraction = fraction + day_diff(date1,date2,"ACT")/365;
	  end	  
	end
end

% funkcja 11.2
% in: data1, data2, formatu jak wyzej
% out: ulamek roku (w konwencji ACT/365L) wyznaczony przez te daty

function fraction = act_365l_frac(date1, date2)  
	tm1 = strptime(date1, "%d-%b-%Y");
	tm2 = strptime(date2, "%d-%b-%Y");
	R1=tm1.year+1900; R2=tm2.year+1900;
	
	if months == 12
	  if ((is_leap_year(R1) && (tm1.mon==0 || (tm1.mon==1 && tm1.mday!=29))) || (is_leap_year(R2) && (tm2.mon>1 || (tm2.mon==1 && tm1.mday==29))))
	  fraction = day_diff(date1,date2,"ACT")/366;
	  else
	  fraction = day_diff(date1,date2,"ACT")/365;	  
	  end
	else
	fraction = day_diff(date1,date2,"ACT")/(365+is_leap_year(R2));
	end
end

% funkcja 11.3
% in: date1, date2, formatu jak wyzej (poczatek i koniec okresu, dla ktorego chcemy poznac date roku)
% date_start1, date_end1 (poczatek i koniec okresu odsetkowego, w ktorym lezy data1), formatu jak wyzej, 
% date_start2, date_end2 (poczatek i koniec okresu odsetkowego, w ktorym lezy data2), formatu jak wyzej,
% period1 period2, period_mid (dlugosci odpowiednio: okresu w ktorym lezy data1, okresu w ktorym lezy data2 i okresu pomiedzy tymi dwoma okresami,
% 	wyrazone w miesiacach; dla period_mid nalezy podac 0 jesli okresy, w ktorych lezy data1 i data 2 sie pokrywaja lub ze soba bezsposrednio sasiaduja)
% out: ulamek roku (w konwencji ACT/ACT ICMA) wyznaczony przez te daty
%
% ta funkcja nie jest dostepna przez year_frac(), a tym samym przez payments(), gdyz potrzebuje dodatkowych parametrow

function fraction = act_act_icma_frac(date1, date2, date_start1, date_end1, date_start2, date_end2, period1, period2, period_mid) 
	if (day_diff(date_start1,date_start2,"ACT")==0) %przypadek, gdy poczatek i koniec leza w tym samym okresie odsetkowym
	  fraction = day_diff(date1,date2,"ACT")/(12/period1*day_diff(date_start1,date_end1,"ACT"))
	else
	  fraction = day_diff(date1,date_end1,"ACT")/(12/period1*day_diff(date_start1,date_end1,"ACT"))+period_mid/12+day_diff(date_start2,date2,"ACT")/(12/period2*day_diff(date_start2,date_end2,"ACT"));
	end
end

% funkcja 11.4
% in: data1, data2, maturity date formatu jak wyzej
% out: ulamek roku (w konwencji 30E/360 ISDA) wyznaczony przez te daty
%
% ta funkcja nie jest dostepna przez year_frac(), a tym samym przez payments(), gdyz potrzebuje dodatkowego parametru

function fraction = isda_30e_360_frac(date1, date2, maturity)
	tm1 = strptime(date1, "%d-%b-%Y");
	tm2 = strptime(date2, "%d-%b-%Y");
	R1=tm1.year+1900; R2=tm2.year+1900;
	
	if (tm1.mday==eomday(R1,tm1.mon+1)) tm1.mday=30; end
	if (tm2.mday==eomday(R2,tm2.mon+1) && !(strcmp(date2,maturity) && tm2.mon==1)) tm2.mday=30; end	
	fraction = ((tm2.year-tm1.year)*360+(tm2.mon-tm1.mon)*30+tm2.mday-tm1.mday)/360;
end

% funkcja 12 (zmieniona)
% in: data1, data2, formatu jak wyzej, Day Count Convention,
% liczba miesiecy pomiedzy data1 i data2 (lub 0, jesli platnosc jest okreslona w sposob inny niz "za n miesiecy")
% out: ulamek roku wyznaczony przez te daty
% 
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function fraction = year_frac(date1, date2, DCC) 
	if strcmp(DCC,"ACT/365") fraction = day_diff(date1,date2,"ACT")/365;
	elseif strcmp(DCC,"ACT/360") fraction = day_diff(date1,date2,"ACT")/360;
	elseif strcmp(DCC,"ACT/ACT") fraction = act_act_frac(date1,date2);
	elseif strcmp(DCC,"ACT/365L") fraction = act_365l_frac(date1,date2);
	elseif strcmp(DCC,"ACT/ACT AFB") fraction = act_act_afb_frac(date1,date2);
	elseif strcmp(DCC,"30/360") fraction = day_diff(date1,date2,"30")/360;
	elseif strcmp(DCC,"30E/360") fraction = day_diff(date1,date2,"30E")/360;
	elseif strcmp(DCC,"30E+/360") fraction = day_diff(date1,date2,"30E+")/360;

	% ewentualne dalsze konwencje
	end
end

% funkcja 13 (zmieniona)
% in: data1, wektor liczby dni w okresach, cell array identyfikatorow FC, Day Count Convention, konwencja Business Day Adjustment, 
% Contract Payment Offset (liczba dni roboczych od contract date do value date)
% out: lista dat (dni robocze) odpowiadajacych datom wyznaczonym z danych wejsciowych, lista ulamków roku pomiedzy nimi
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual
%
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function [data_list,frac_list] = payments(date1,dwect,FCL,DCC,BDA,CPO)
	n = length(dwect);
	d = cumsum(dwect);
	data_list = cell(n+1,1);
	frac_list = cell(n,1);
	if (CPO==0)		% względy czasu wykonywania
		data_list{1} = date_rolling(date1,FCL,BDA);
		tm1 = strptime(date1, "%d-%b-%Y");
		for i=2:n+1
			data_list{i} = date_rolling(strftime("%d-%b-%Y",localtime(mktime(tm1)+d(i-1)*86400+7200)),FCL,BDA);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	else
		data_list{1} = day_shift2(date_rolling(date1,FCL,BDA),FCL,CPO);
		tm1 = strptime(date1, "%d-%b-%Y");
		for i=2:n+1
			data_list{i} = day_shift2(date_rolling(strftime("%d-%b-%Y",localtime(mktime(tm1)+d(i-1)*86400+7200)),FCL,BDA),FCL,CPO);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	end
end

% funkcja 14 (modyfikacja 13) (zmieniona)
% in: lista dat, cell array identyfikatorow FC, Day Count Convention, konwencja Business Day Adjustment, Contract Payment Offset (liczba dni roboczych
% od contract date do value date)
% out: lista dat (dni robocze) odpowiadajacych datom wzynaczonym z danych wejsciowych, lista ulamków roku pomiedzy nimi
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual
%
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function [data_list,frac_list] = payments2(dates,FCL,DCC,BDA,CPO)
	n = length(dates)-1;
	data_list = cell(n+1,1);
	frac_list = cell(n,1);
	if (CPO==0)		% względy czasu wykonywania
		data_list{1} = date_rolling(dates{1},FCL,BDA);
		for i=2:n+1
			data_list{i} = date_rolling(dates{i},FCL,BDA);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	else
		data_list{1} = day_shift2(date_rolling(dates{1},FCL,BDA),FCL,CPO);
		for i=2:n+1
			data_list{i} = day_shift2(date_rolling(dates{i},FCL,BDA),FCL,CPO);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	
	end
end

% funkcja 15 (modyfikacja 13) (zmieniona)
% in: data1, liczba okresow, liczba dni w okresie, cell array identyfikatorow FC, Day Count Convention, konwencja BDA, 
% Contract Payment Offset (liczba dni roboczych od contract date do value date)
% out: lista dat (dni robocze) odpowiadajacych datom wyznaczonym z danych wejsciowych, lista ulamkow roku pomiedzy nimi
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual
%
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function [data_list,frac_list] = payments3(date1,n,d,FCL,DCC,BDA,CPO)
	[data_list,frac_list] = payments(date1,repmat(d,1,n),FCL,DCC,BDA,CPO);
end

% funkcja 16 (modyfikacja 13) (zmieniona)
% in: data1, wektor liczby miesiecy w okresach, cell array identyfikatorow FC, Day Count Convention, 
% konwencja Business Day Adjustment, Contract Payment Offset (liczba dni roboczych od contract date do value date), 
% parametr End of Month Adjustment +/-1 ("w przod"/"w tyl") do add_months()
% out: lista dat (dni robocze) odpowiadajacych datom wyznaczonym z danych wejsciowych, lista ulamkow roku pomiedzy nimi
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual
%
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function [data_list,frac_list] = payments4(date1,mwect,FCL,DCC,BDA,CPO,EMA)
	n = length(mwect);
	d = cumsum(mwect);
	data_list = cell(n+1,1);
	frac_list = cell(n,1);
	ad_mon = add_months(date1,mwect,EMA);
	if (CPO==0)		% wzgledy czasu wykonywania
		data_list{1} = date_rolling(date1,FCL,BDA);
		for i=2:n+1
			data_list{i} = date_rolling(ad_mon{i},FCL,BDA);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	else
		data_list{1} = day_shift2(date_rolling(date1,FCL,BDA),FCL,CPO);
		for i=2:n+1
			data_list{i} = day_shift2(date_rolling(ad_mon{i},FCL,BDA),FCL,CPO);
			frac_list{i-1} = year_frac(data_list{i-1},data_list{i},DCC);
		end	
	end
end

% funkcja 17 (modyfikacja 13) (zmieniona)
% in: data1, liczba okresow, liczba miesiecy w okresie, cell array identyfikatorow FC, Day Count Convention, 
% konwencja Business Day Adjustment, Contract Payment Offset (liczba dni roboczych od contract date do value date), 
% parametr End of Month Adjustment +/-1 ("w przod"/"w tyl") do add_months()
% out: lista dat (dni robocze) odpowiadajacych datom wyznaczonym z danych wejsciowych, lista ulamkow roku pomiedzy nimi
%
% identyfikatory BDA:
% sfbd = Standard Following Business Day, mfbd = Modified Following Business Day
% spbd = Standard Previous Business Day, mpbd = Modified Previous Business Day
% seom = Standard End of Month, actu = Actual
%
% identyfikatory DCC: ACT/365, ACT/360, ACT/ACT, ACT/365L, ACT/ACT AFB, 30/360, 30E/360, 30E+/360   

function [data_list,frac_list] = payments5(date1,n,m,FCL,DCC,BDA,CPO, EMA)
	[data_list,frac_list] = payments4(date1,repmat(m,1,n),FCL,DCC,BDA,CPO, EMA);
end
