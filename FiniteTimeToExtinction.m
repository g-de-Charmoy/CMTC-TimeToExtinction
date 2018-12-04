%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Approximation of Finite Time to Extinction for CTMC Ebola Outbreak Model using 10000 Sample Paths
%By Gabriel de Charmoy
%Hons Student (Applied Mathematics) - University of KwaZulu-Natal, RSA
%13 October 2018



clear all
close all
clc;
set(0,'DefaultAxesFontSize', 18)%Increases axes labels
set(gca,'fontsize',18);%Changes font size in a figure.
%-------------------------------------------------------------------------
%Parameter values of the model
time=0.6; %time
p=0.1;  %Quarantine factor E
q=0.8;  %Quartine factor I
r=0.5;  %Recover factor
k=0.025; %infection decay factor %0.0275
tau=9999; %decay delay
beta =1.225e-4; %2e-4
nd = 0.05;
sigma = (1-p)/11;
phi = r/4;
gamma = (1-q)/5.3;
alphaQ = (1-r)/4;
alphaD = 0.25;
thetaI = q/5.3;
thetaE = p/11;
N=10000;
%-------------------------------------------------------------------------
%Initial Conditions
init = [N, 0, 1, 0, 0, 0, 0];
%-------------------------------------------------------------------------
%Three sample paths for CTMC (Gillespie's algorithm)
z=0;
ext=10000;
fatality=zeros(ext,1);
deaths=zeros(ext,1);
cases=zeros(ext,1);
maxInfectious=zeros(ext,1);
end_times=zeros(ext,1);  % array to store ending times

for s=1:ext
i=0;
c=0;
d=0;
clear t S E I Q R D
j=1;
S=init(1); E=init(2); I=init(3);Q=init(4); R=init(5); D=init(6); 
t=0; %Starting at zero
%-------------------------------------------------------------------------
while   (E+I+D)>0 && t<time    % Stop if it hits zero or at end time  
     
randm1=rand; %Uniform random number
randm2=rand; %Uniform random number
 %-------------------------------------------------------------------------  
    

 %-------------------------------------------------------------------------    
 %Defining all possible events at time (t+1)
   ev1 = beta*exp(-k*t)*S*(I+nd*D);
    ev2 = ev1 + sigma*E; 
    ev3 = ev2 + thetaE*E;
    ev4 = ev3 + thetaI*I;
    ev5 = ev4 + phi*Q;
    ev6 = ev5 + gamma*I;
    ev7 = ev6 + alphaD*D;
    ev8 = ev7 + alphaQ*Q;
 %-------------------------------------------------------------------------  
   t=t+abs((log(0.8+randm1*0.1)))/ev8; % Time to next event
   randm2=randm2*ev8;
 %-------------------------------------------------------------------------   
 
 
%Computing probabilities of each event occuring at time (t+1)
     if randm2<ev1;  
        S = S-1;
        E = E+1;
        I = I;
        Q = Q;
        R = R;
        D = D;
        c=c+1;
    elseif randm2<ev2
        E = E-1;
        I = I+1;
    elseif  randm2<ev3
        E = E-1;
        Q = Q+1;
    elseif  randm2<ev4
        I = I-1;
        Q = Q+1;
   elseif  randm2<ev5
        Q = Q-1;
        R = R+1;
   elseif  randm2<ev6
        I = I-1;
        D = D+1;
    elseif  randm2<ev7
        D = D-1;
        d=d+1;
     else  
        Q = Q-1;
        d=d+1;
     end
    j=j+1;
    
    i=max(i,I+D+E);
end
if (I+E+D)==0
    z=z+1;
else
    z=z;
end
deaths(s)=d;
cases(s)=c;

maxInfectious(s)=i;
end_times(s)=t;
end

z;  %number of sample paths that hit zero out of 8000 sample paths
probext=z/ext  %approximate probability of extinction
%hist(maxInfectious,100)
%hist(maxInfectious,100)
%xtickformat('%.1f')
hist(end_times,100)
title('Peak Size of Epidemic (I+E+D)', 'fontsize', 18)
%xlabel('Extinction Time (Days)')
xlabel('Number of infectives ')
%xlabel('Number of Cases')
ylabel('Probability (x10-4)')


%xlim([0 max(maxInfectious)]) % Max Size of infection
%xlim([0 max(fatality)])%case count
%xlim([0 time]) %Extinction Time
%ylim([0 1500]) %y-axis for all three present
%ylim([0 4000]) %y-axis for X12
end_times; % Shows end times per simulation
T=mean(end_times) %To calculate the mean extinction time...
Average_Peak_Size_of_Infection=mean(maxInfectious)
Average_number_of_Cases=mean(cases)
Average_falaity_rate=mean(deaths)/mean(cases)
