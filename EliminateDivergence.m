 clc;clear;

% y1 = textread('1 Voltage update1.dat','','headerlines',0); %Read file:Time	Uccp	Usource	Uc	Qc	I_ec
%  y1 = textread('Grid0D       21004 .dat','','headerlines',0); %Read file: Time	I	Psource	Pccp	Pc
 y1 = textread('1 EC update1.dat','','headerlines',0); %Read file: Time	Uccp	Usource	Uc	Qc	I_ec  Psource	Pccp  Pc    Pfield

t_s = 5.0e-11; %????
dt = t_s;
%  t_s = 1/60.0e+6; %????
t_start = 0.0e-6+t_s; %
t_end = 1.0e-6;     %????
t = t_start : t_s : t_end;
num = length(t);
nn  = 10;
num_n = floor(num/nn);
% y = 1.5*sin(2*pi*5*t)+3*sin(2*pi*20*t)+randn(1,length(t));  %????
% for 1:num
%     y =  y1(:,2:num);
% end 
T = 1/60.0E6;
nstart = floor(t_start/t_s);
t_temp = y1(nstart:(nstart+num-1),1);
yQ = y1(nstart:(nstart+num-1),5);
yUccp = y1(nstart:(nstart+num-1),2);
yIec = y1(nstart:(nstart+num-1),6);
yPccp = y1(nstart:(nstart+num-1),8);
yPc = y1(nstart:(nstart+num-1),9);
yPfield = y1(nstart:(nstart+num-1),10);
 yEfield = y1(nstart:(nstart+num-1),11);

for n = 2:num_n-1
    if yUccp(n)<1
%         yUccp(n) = yUccp(n-1);
          yUccp(n) =1;
     end
    tfigure(n) =  t_temp(nn*n) ;
%     R(n) =  ( yPccpave(n) - yPfieldave(n))/(yIave(n)*yIave(n));
    C(n) =   2*yEfield(n)/(yUccp(n)*yUccp(n));
end

plot(tfigure,C);