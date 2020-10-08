 clc;clear;

% y1 = textread('1 Voltage update1.dat','','headerlines',0); %Read file:Time	Uccp	Usource	Uc	Qc	I_ec
%  y1 = textread('Grid0D       21004 .dat','','headerlines',0); %Read file: Time	I	Psource	Pccp	Pc
 y1 = textread('1 EC update1.dat','','headerlines',0); %Read file: Time	Uccp	Usource	Uc	Qc	I_ec  Psource	Pccp  Pc    Pfield
%Ef     \Phi	P_heat	P_coll	P_boundary	P_heat-P_coll

t_s = 5.0e-11; %????
dt = t_s;
%  t_s = 1/60.0e+6; %????
t_start =1.62e-6+t_s; %
t_end =1.66e-6;     %????
t = t_start : t_s : t_end;
num = length(t);
nn  = 1;
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
 
 Tnum = 1;
for n =1:num_n-1
%     if n==1
%         E_0(n)= yEfield(n);
%     else
%     
%      E_0(n)=E_0(n-1);
         
    tfigure(n) =  t_temp(nn*n);
    if (yUccp(n)*yUccp(n+1)<0)
          if(abs(yUccp(n))>abs(yUccp(n+1)))
              Tnum= n+1;
          else
              Tnum= n;
          end

    end
          Efield0(n)  = yEfield(n);
          Efield_0(n) = yEfield(Tnum);
          Efield(n) = yEfield(n)- Efield_0(n);
          Pccp(n) = yPccp(n);
          P_Efield(n) = yPfield(n);
          U_ccp(n) =  yUccp(n);
          U_cccp(n) =  U_ccp(n) * P_Efield(n)/Pccp(n);
      
   C1(n) =   2*yEfield(n) /U_ccp(n)/U_ccp(n) ;
   C2(n) =   2*Efield(n) /U_ccp(n)/U_ccp(n) ;
   C3(n) =   2*Efield(n) /U_cccp(n)/U_cccp(n);
%    C4(n) =   2*Efield_0(n) /U_cccp(n)/U_cccp(n);
end

% Efield(n)
% figure(3,1);
% plot(tfigure, U_ccp,'R', tfigure, U_ccp_c,'B');
%  plot(tfigure, Efield,'R', tfigure,Efield0,'B')%, tfigure,U_ccp,'G');
plot(tfigure, C1,'R', tfigure,C2,'B',tfigure,C3,'k')%, tfigure,U_ccp,'G');
%   plot(tfigure,C2,'B')%, tfigure,U_ccp,'G');
% figure(3,2);
% plot(tfigure, E_field);
% figure(3,3);
% hold on;
%  plot(tfigure,C);
