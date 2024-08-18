clc
clear all 
close all


%% Load flow Variables
basekv=12.66;
basemva = 100;
tic
base_impedance=(basekv)^2/basemva; %base impedance=(basekv)^2/basemva
accuracy = 0.001;
maxiter = 1000;
x=.001;
s=1/base_impedance;

% Hourly load (in % w.r.t 24 hours)
loadddlevel2=[65 60 56 57 58 58.5 59 62 67 69.5 72 79 88 93 98 100 97 93 88 85.5 82 73 62 58];


%        IEEE 33-BUS TEST SYSTEM (American Electric Power)
%        Bus Bus  Voltage Angle   ---Load---- -------Generator----- Static Mvar
%        No  code Mag.    Degree  MW    Mvar  MW  Mvar Qmin Qmax    +Qc/-Ql
busdata=[1   1    1    0      0   0      0    0    0   0    0
         2   0    1    0       100*x   60*x       0    0    0   0    0
         3   0    1   0        90*x    40*x      0    0    0    0       0
         4   0    1   0        120*x   80*x      0    0    0    0       0
         5   0    1    0        60*x    30*x      0    0    0    0       0
         6   0    1    0        60*x    20*x      0    0    0    0       0
         7   0    1    0        200*x   100*x     0    0    0    0       0       
         8   0    1    0        200*x   100*x     0    0    0    0       0  
         9   0    1    0        60*x    20*x      0    0    0    0       0
         10   0    1    0        60*x    20*x      0    0    0    0       0
         11   0    1    0        45*x    30*x      0    0    0    0       0
         12   0    1    0        60*x    35*x      0    0    0    0       0
         13   0    1    0        60*x    35*x      0    0    0    0       0
         14   0    1    0        120*x   80*x      0    0    0    0       0
         15   0    1    0        60*x    10*x      0    0    0    0       0
         16   0    1    0        60*x    20*x      0    0    0    0       0
         17   0    1    0        60*x    20*x      0    0    0    0       0
         18   0    1    0        90*x    40*x      0    0    0    0       0
         19   0    1    0        90*x    40*x      0    0    0    0       0
         20   0    1    0        90*x    40*x      0    0    0    0       0 
         21   0    1    0        90*x    40*x      0    0    0    0       0
         22   0    1    0        90*x    40*x      0    0    0    0       0
         23   0    1    0        90*x    50*x      0    0    0    0       0
         24   0    1    0        420*x   200*x     0    0    0    0       0
         25   0    1    0        420*x   200*x     0   0    0    0       0
         26   0    1    0        60*x    25*x      0    0    0    0       0
         27   0    1    0        60*x    25*x      0    0    0    0       0
         28   0    1    0        60*x    20*x      0    0    0    0       0
         29   0    1    0        120*x   70*x      0    0    0    0       0
         30   0    1    0        200*x   600*x     0     0    0    0       0
         31   0    1    0        150*x   70*x      0    0    0    0       0  
         32   0    1    0        210*x   100*x     0    0    0    0       0
         33   0    1    0        60*x    40*x      0    0    0    0       0];

linedata=[1	2	0.0922*s 	0.047*s   0   1
2	3	0.493*s 	0.2511*s   0   1
3	4	0.366*s 	0.1864*s   0   1
4	5	0.3811*s 	0.1941*s   0   1
5	6	0.819*s 	0.707*s   0   1
6	7	0.1872*s 	0.6188*s   0   1
7	8	0.7114*s 	0.2351*s   0   1
8	9	1.03*s 	    0.74*s   0   1
9	10	1.044*s 	0.74*s   0   1
10	11	0.1966*s 	0.065*s   0   1
11	12	0.3744*s 	0.1238*s   0   1
12	13	1.468*s 	1.155*s   0   1
13	14	0.5416*s 	0.7129*s   0   1
14	15	0.591*s 	0.526*s   0   1
15	16	0.7463*s 	0.545*s   0   1
16	17	1.289*s 	1.721*s   0   1
17	18	0.732*s 	0.574*s   0   1
2	19	0.164*s 	0.1565*s   0   1
19	20	1.5042*s 	1.3554*s   0   1
20	21	0.4095*s 	0.4784*s   0   1
21	22	0.7089*s 	0.9373*s   0   1
3	23	0.4512*s 	0.3083*s   0   1
23	24	0.898*s 	0.7091*s   0   1
24	25	0.896*s 	0.7011*s   0   1
6	26	0.203*s 	0.1034*s   0   1
26	27	0.2842*s 	0.1447*s   0   1
27	28	1.059*s 	0.9337*s   0   1
28	29	0.8042*s 	0.7006*s   0   1
29	30	0.5075*s 	0.2585*s   0   1
30	31	0.9744*s 	0.963*s   0   1
31	32	0.3105*s 	0.3619*s   0   1
32	33	0.341*s 	0.5302*s   0   1];


%% First three are higher limits of Constant generators with the las three being positional indexes
highlimit =[1.48 1.48 1.48 33 33 33];
lowlimit=[0.1 .1 0.1 2 2 2];

dim=6;
Run_no=1;         % Number of independent runs 
Particles_no=50;   % Number of particles
Max_iterationsss=100; % Maximum number of iterations


for irun=1:Run_no

Ceq1=zeros(1,dim);   Ceq1_fit=inf; 
Ceq2=zeros(1,dim);   Ceq2_fit=inf; 
Ceq3=zeros(1,dim);   Ceq3_fit=inf; 
Ceq4=zeros(1,dim);   Ceq4_fit=inf;

%C=initialization(Particles_no,dim,ub,lb);

for hhjk=1:50
 Cpop (hhjk,:) = [(1.48)*rand(1,3) randi([2,33],1,3)];

end



Itesdr=0; Veo=1;

as1=2;
as2=1;
GsP=0.5;
dsx(1)=.7;


a=0.7;
for i=1:500
    dsx(i+1)=sin((a*pi)/dsx(i));
    G(i)=((dsx(i)+1)*100)/2;
end
%normalize it from [-1 1] to [0 1]
a=-1; b=1; c=0; d=1;
dsx=((dsx-a)*(d-c))/(b-a);

%for i=1:Max_iter
 %    dsx(i+1) = 2.3*dsx(i)^2*sin(pi*dsx(i));
  %   G(i)=(dsx(i))*100;
 %end

while Itesdr<Max_iterationsss
   
        
        
          
        for ikodh=1:50
     for jmm=1:dim
    
if Cpop(ikodh,jmm)>highlimit(1,jmm) || Cpop(ikodh,jmm)<lowlimit(1,jmm)
     Cpop(ikodh,:)=[(1.48)*rand(1,3) randi([2,33],1,3) ]; %randi([2,33],1,2)
                        
               
          else Cpop(ikodh,jmm)= Cpop(ikodh,jmm);
end
     end
 end  
    
 for ikodh=1:50
   
     busdata= [1   1    1    0      0   0      0    0    0   0    0
         2   0    1    0       100*x   60*x       0    0    0   0    0
         3   0    1   0        90*x    40*x      0    0    0    0       0
         4   0    1   0        120*x   80*x      0    0    0    0       0
         5   0    1    0        60*x    30*x      0    0    0    0       0
         6   0    1    0        60*x    20*x      0    0    0    0       0
         7   0    1    0        200*x   100*x     0    0    0    0       0       
         8   0    1    0        200*x   100*x     0    0    0    0       0  
         9   0    1    0        60*x    20*x      0    0    0    0       0
         10   0    1    0        60*x    20*x      0    0    0    0       0
         11   0    1    0        45*x    30*x      0    0    0    0       0
         12   0    1    0        60*x    35*x      0    0    0    0       0
         13   0    1    0        60*x    35*x      0    0    0    0       0
         14   0    1    0        120*x   80*x      0    0    0    0       0
         15   0    1    0        60*x    10*x      0    0    0    0       0
         16   0    1    0        60*x    20*x      0    0    0    0       0
         17   0    1    0        60*x    20*x      0    0    0    0       0
         18   0    1    0        90*x    40*x      0    0    0    0       0
         19   0    1    0        90*x    40*x      0    0    0    0       0
         20   0    1    0        90*x    40*x      0    0    0    0       0 
         21   0    1    0        90*x    40*x      0    0    0    0       0
         22   0    1    0        90*x    40*x      0    0    0    0       0
         23   0    1    0        90*x    50*x      0    0    0    0       0
         24   0    1    0        420*x   200*x     0    0    0    0       0
         25   0    1    0        420*x   200*x     0   0    0    0       0
         26   0    1    0        60*x    25*x      0    0    0    0       0
         27   0    1    0        60*x    25*x      0    0    0    0       0
         28   0    1    0        60*x    20*x      0    0    0    0       0
         29   0    1    0        120*x   70*x      0    0    0    0       0
         30   0    1    0        200*x   600*x     0     0    0    0       0
         31   0    1    0        150*x   70*x      0    0    0    0       0  
         32   0    1    0        210*x   100*x     0    0    0    0       0
         33   0    1    0        60*x    40*x      0    0    0    0       0];
     
    
    
    busdata(ceil(Cpop(ikodh,4)),7)=.9*Cpop(ikodh,1);
    busdata(ceil(Cpop(ikodh,4)),8)=((Cpop(ikodh,1))^2 - (0.9*Cpop(ikodh,1))^2)^0.5;
    busdata(ceil(Cpop(ikodh,5)),7)=.9*Cpop(ikodh,2);
    busdata(ceil(Cpop(ikodh,5)),8)=((Cpop(ikodh,2))^2 - (0.9*Cpop(ikodh,2))^2)^0.5;
    busdata(ceil(Cpop(ikodh,6)),7)=.9*Cpop(ikodh,3);
    busdata(ceil(Cpop(ikodh,6)),8)=((Cpop(ikodh,3))^2 - (0.9*Cpop(ikodh,3))^2)^0.5;
    %busdata(ceil(Cpop(ikodh,6)),7)=Cpop(ikodh,2);
    %busdata(ceil(Cpop(ikodh,7)),7)=Cpop(ikodh,3);
    %busdata(ceil(Cpop(ikodh,8)),7)=Cpop(ikodh,4);
    
    busdata(ceil(Cpop(ikodh,5)),2)=0;
    busdata(ceil(Cpop(ikodh,4)),2)=0;
    busdata(ceil(Cpop(ikodh,6)),2)=0;
    %busdata(ceil(Cpop(ikodh,8)),2)=2;
lfybusold;
lfnewton12;
lineflow2;

%% Constraints
deviation= sum(abs(ones(1,33)-Vm));
if min(Vm)<.95|| max(Vm)>1.05
     fitness(ikodh,:)=1000;
     
    
else fitness(ikodh,:)=((Loss+deviation));
 end
 end    
          
        %Flag4ub=C(i,:)>ub;
        %Flag4lb=C(i,:)<lb;
        %C(i,:)=(C(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;         
          
       % fitness(i)=fobj(C(i,:));
      for gi=1:size(Cpop,1)
        if fitness(gi)<Ceq1_fit 
              Ceq1_fit=fitness(gi);  Ceq1=Cpop(gi,:);
        elseif fitness(gi)>Ceq1_fit && fitness(gi)<Ceq2_fit  
              Ceq2_fit=fitness(gi);  Ceq2=Cpop(gi,:);              
        elseif fitness(gi)>Ceq1_fit && fitness(gi)>Ceq2_fit && fitness(gi)<Ceq3_fit
              Ceq3_fit=fitness(gi);  Ceq3=Cpop(gi,:);
        elseif fitness(gi)>Ceq1_fit && fitness(gi)>Ceq2_fit && fitness(gi)>Ceq3_fit && fitness(gi)<Ceq4_fit
              Ceq4_fit=fitness(gi);  Ceq4=Cpop(gi,:);
                         
        end
      end
      
%---------------- Memory saving-------------------   
      if Itesdr==0
        fit_old=fitness;  C_old=Cpop;
      end
    
     for gi=1:Particles_no
         if fit_old(gi)<fitness(gi)
             fitness(gi)=fit_old(gi); Cpop(gi,:)=C_old(gi,:);
         end
     end

    C_old=Cpop;  fit_old=fitness;
%-------------------------------------------------
       
Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;                              % averaged candidate 
C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];                     % Equilibrium pool

 
 t=(1-Itesdr/Max_iterationsss)^(as2*Itesdr/Max_iterationsss);                      % Eq (9)

 
    for gi=1:Particles_no
           lambda=rand(1,dim);                                % lambda in Eq(11)
           r=rand(1,dim);                                     % r in Eq(11)  
           Ceq=C_pool(randi(size(C_pool,1)),:);               % random selection of one candidate from the pool
           Feo=abs(as1*sign(r-0.5).*(exp(-lambda.*t)-1));             % Eq(11)
           r1=rand();
           
           r2=rand();                              % r1 and r2 in Eq(15)
           GCP=0.5*r1*ones(1,dim)*(r2>=GsP);                   % Eq(15)
           Geo0=GCP.*(Ceq-lambda.*Cpop(gi,:));                      % Eq(14)
           Geo=Geo0.*Feo;                                           % Eq(13)
           Cpop(gi,:)=Ceq+(Cpop(gi,:)-Ceq).*Feo+(Geo./lambda.*Veo).*(1-Feo);   % Eq(16)                                                             
    end
 
       Itesdr=Itesdr+1;  
       Convergence_curve(Itesdr)=Ceq1_fit; 
       Ceqfit_run(irun)=Ceq1_fit;
fprintf('eoA|%5.0f -----> %9.16f\n',Itesdr,Ceq1_fit)

end


display(['Run no : ', num2str(irun)]);
display(['The best solution obtained by CEO is : ', num2str(Ceq1,10)]);
display(['The best optimal value of the objective funciton found by EO is : ', num2str(Ceq1_fit,10)]);
disp(sprintf('--------------------------------------'));
end


Ave=mean(Ceqfit_run);
Sd=std(Ceqfit_run);

 busdata=[1   1    1    0      0   0      0    0    0   0    0
         2   0    1    0       100*x   60*x       0    0    0   0    0
         3   0    1   0        90*x    40*x      0    0    0    0       0
         4   0    1   0        120*x   80*x      0    0    0    0       0
         5   0    1    0        60*x    30*x      0    0    0    0       0
         6   0    1    0        60*x    20*x      0    0    0    0       0
         7   0    1    0        200*x   100*x     0    0    0    0       0       
         8   0    1    0        200*x   100*x     0    0    0    0       0  
         9   0    1    0        60*x    20*x      0    0    0    0       0
         10   0    1    0        60*x    20*x      0    0    0    0       0
         11   0    1    0        45*x    30*x      0    0    0    0       0
         12   0    1    0        60*x    35*x      0    0    0    0       0
         13   0    1    0        60*x    35*x      0    0    0    0       0
         14   0    1    0        120*x   80*x      0    0    0    0       0
         15   0    1    0        60*x    10*x      0    0    0    0       0
         16   0    1    0        60*x    20*x      0    0    0    0       0
         17   0    1    0        60*x    20*x      0    0    0    0       0
         18   0    1    0        90*x    40*x      0    0    0    0       0
         19   0    1    0        90*x    40*x      0    0    0    0       0
         20   0    1    0        90*x    40*x      0    0    0    0       0 
         21   0    1    0        90*x    40*x      0    0    0    0       0
         22   0    1    0        90*x    40*x      0    0    0    0       0
         23   0    1    0        90*x    50*x      0    0    0    0       0
         24   0    1    0        420*x   200*x     0    0    0    0       0
         25   0    1    0        420*x   200*x     0   0    0    0       0
         26   0    1    0        60*x    25*x      0    0    0    0       0
         27   0    1    0        60*x    25*x      0    0    0    0       0
         28   0    1    0        60*x    20*x      0    0    0    0       0
         29   0    1    0        120*x   70*x      0    0    0    0       0
         30   0    1    0        200*x   600*x     0     0    0    0       0
         31   0    1    0        150*x   70*x      0    0    0    0       0  
         32   0    1    0        210*x   100*x     0    0    0    0       0
         33   0    1    0        60*x    40*x      0    0    0    0       0];

linedata=[1	2	0.0922*s 	0.047*s   0   1
2	3	0.493*s 	0.2511*s   0   1
3	4	0.366*s 	0.1864*s   0   1
4	5	0.3811*s 	0.1941*s   0   1
5	6	0.819*s 	0.707*s   0   1
6	7	0.1872*s 	0.6188*s   0   1
7	8	0.7114*s 	0.2351*s   0   1
8	9	1.03*s 	    0.74*s   0   1
9	10	1.044*s 	0.74*s   0   1
10	11	0.1966*s 	0.065*s   0   1
11	12	0.3744*s 	0.1238*s   0   1
12	13	1.468*s 	1.155*s   0   1
13	14	0.5416*s 	0.7129*s   0   1
14	15	0.591*s 	0.526*s   0   1
15	16	0.7463*s 	0.545*s   0   1
16	17	1.289*s 	1.721*s   0   1
17	18	0.732*s 	0.574*s   0   1
2	19	0.164*s 	0.1565*s   0   1
19	20	1.5042*s 	1.3554*s   0   1
20	21	0.4095*s 	0.4784*s   0   1
21	22	0.7089*s 	0.9373*s   0   1
3	23	0.4512*s 	0.3083*s   0   1
23	24	0.898*s 	0.7091*s   0   1
24	25	0.896*s 	0.7011*s   0   1
6	26	0.203*s 	0.1034*s   0   1
26	27	0.2842*s 	0.1447*s   0   1
27	28	1.059*s 	0.9337*s   0   1
28	29	0.8042*s 	0.7006*s   0   1
29	30	0.5075*s 	0.2585*s   0   1
30	31	0.9744*s 	0.963*s   0   1
31	32	0.3105*s 	0.3619*s   0   1
32	33	0.341*s 	0.5302*s   0   1];


busdata(ceil(Ceq1(1,4)),7)=0.9*Ceq1(1,1);
   busdata(ceil(Ceq1(1,4)),8)=(   (Ceq1(1,1))^2 - (0.9*Ceq1(1,1))^2)^0.5;
  
 busdata(ceil(Ceq1(1,5)),7)=0.9*Ceq1(1,2);
   busdata(ceil(Ceq1(1,5)),8)=(   (Ceq1(1,2))^2 - (0.9*Ceq1(1,2))^2)^0.5;
 busdata(ceil(Ceq1(1,6)),7)=0.9*Ceq1(1,3);
   busdata(ceil(Ceq1(1,6)),8)=(   (Ceq1(1,3))^2 - (0.9*Ceq1(1,3))^2)^0.5;
   lfybusold
    lfnewton12
    busout
    lineflow2
    Loss
    deviation= sum(abs(ones(1,33)-Vm));
et=toc
 hh=[Loss deviation Ceq1 Ceq1_fit et]'

