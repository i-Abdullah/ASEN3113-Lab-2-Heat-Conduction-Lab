%% housekeeping

clear
clc
close all


%% read data

% add data files path
addpath('./Data');

% read data
Steel = importdata('Steel_20V_185mA.txt'); % the steel sepcimen data
Brass = importdata('Brass_30V_285mA.txt'); % the brass sepcimen data
Alum = importdata('Aluminum_25V_240mA.txt'); % the aluminum specimen data

% voltage
V_Steel = 20;
V_Brass = 30;
V_Alum = 25;

% current
A_Steel = 185;
A_Brass = 285;
A_Alum = 240;


% k : tehrmaol conductivity:

k_steel = 16.2; % W / m.k;
k_Alum = 130;
k_Brass = 115;

% convert to inches:

k_steel = k_steel*0.0254; % W / m.k;
k_Alum = k_Alum*0.0254;
k_Brass = k_Brass*0.0254;



% parse the data:


% temp is in C?

time_Steel = Steel.data(:,1); % time
TC1_Steel = Steel.data(:,2); % temp reading from first thermo couple
TC2_Steel = Steel.data(:,3); % temp reading from 2nd thermo couple
TC3_Steel = Steel.data(:,4); % temp reading from 3rd thermo couple
TC4_Steel = Steel.data(:,5); 
TC5_Steel = Steel.data(:,6); 
TC6_Steel = Steel.data(:,7); 
TC7_Steel = Steel.data(:,8); 
TC8_Steel = Steel.data(:,9);
TC_Steel = Steel.data(:,(2:9)); % all tehrmocouples temp.


time_Brass = Brass.data(:,1); % time
TC1_Brass = Brass.data(:,2); % temp reading from first thermo couple
TC2_Brass = Brass.data(:,3); % temp reading from 2nd thermo couple
TC3_Brass = Brass.data(:,4); % temp reading from 3rd thermo couple
TC4_Brass = Brass.data(:,5); 
TC5_Brass = Brass.data(:,6); 
TC6_Brass = Brass.data(:,7); 
TC7_Brass = Brass.data(:,8); 
TC8_Brass = Brass.data(:,9);
TC_Brass = Brass.data(:,(2:9)); % all tehrmocouples temp.


time_Alum = Alum.data(:,1); % time
TC1_Alum = Alum.data(:,2); % temp reading from first thermo couple
TC2_Alum = Alum.data(:,3); % temp reading from 2nd thermo couple
TC3_Alum = Alum.data(:,4); % temp reading from 3rd thermo couple
TC4_Alum = Alum.data(:,5); 
TC5_Alum = Alum.data(:,6); 
TC6_Alum = Alum.data(:,7); 
TC7_Alum = Alum.data(:,8); 
TC8_Alum = Alum.data(:,9);
TC_Alum = Alum.data(:,(2:9)); % all tehrmocouples temp.

diameter = 1; % in meters 

%% question 1:

% estimate T0:

% T0 is temprature at the cold end, it should be constant through out the
% whole expirement, but we will take it for end time
% to estimate H

% assuming temprature along the rod length changes linearly


% get all temp values for t = end


% we assume cold end starts at x = 0, 0.5 inches, there's a thermocouple,
% another 0.5 inches (i.e. @ x = 1.0)



loc = [ 0.5:0.5:0.5*8 ];

L = loc(end)+1;

endTemp_Steel = TC_Steel(end,:); % temprature of all thermocouples for steel at t = end;
endTemp_Brass = TC_Brass(end,:); % temprature of all thermocouples for steel at t = end;
endTemp_Alum = TC_Alum(end,:); % temprature of all thermocouples for steel at t = end;

[xData, yData] = prepareCurveData( loc, endTemp_Steel );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
T0_Steel = fitresult.p2;
H_Steel_exp = fitresult.p1;

plot( loc, endTemp_Steel , 'sr','MarkerSize',7,'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6] )
hold on
fit_plot = plot(fitresult,'k-');
set(fit_plot, 'LineWidth',2);
hold on

[xData, yData] = prepareCurveData( loc, endTemp_Brass );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
T0_Brass = fitresult.p2;
H_Brass_exp = fitresult.p1;

plot( loc, endTemp_Brass , 'sc','MarkerSize',7,'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[ 0 0.5 1 ] )
hold on
fit_plot = plot(fitresult,'k-');
set(fit_plot, 'LineWidth',2);
hold on



[xData, yData] = prepareCurveData( loc, endTemp_Alum );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
T0_Alum = fitresult.p2;
H_Alum_exp = fitresult.p1;


plot( loc, endTemp_Alum , 's','MarkerSize',7,'MarkerEdgeColor',[1 0.5 0],...
    'MarkerFaceColor',[ 1 0.9 0 ] )
hold on
fit_plot = plot(fitresult,'k-');
set(fit_plot, 'LineWidth',2);
hold on



legend('Steel temp measurements','Best linear fit',...
    'Brass temp measurements','Best linear fit',...
    'Aluminum temp measurements','Best linear fit','Location','NorthWest') 


grid minor


title(['Temperature profile for the steady state experimental data'])
xlabel('Location along the bar [in]')
ylabel('Temp in [C]')


% analytical H, assume Qdot = P = IV

H_Steel_anal = (V_Steel*A_Steel*10^-3) / ( pi*(diameter/2)^2*k_steel ); % analytical H.
H_Alum_anal = (V_Alum*A_Alum*10^-3) / ( pi*(diameter/2)^2*k_Alum ); % analytical H.
H_Brass_anal = (V_Brass*A_Brass*10^-3) / ( pi*(diameter/2)^2*k_Brass ); % analytical H.


H_Steel = [ H_Steel_exp ; H_Steel_anal]./0.0254;
H_Alum = [ H_Alum_exp ; H_Alum_anal ]./0.0254;
H_Brass = [ H_Brass_exp ; H_Brass_anal]./0.0254;
Names = { 'Experimental';'Analytical' };

table(Names,H_Steel,H_Alum,H_Brass)


%% question 2:

%{
Q2 wants you to use fourier series to estimate transient + steady state
sloutions at the same time, thus some fourier series terms must be used.
This sloution will either be analytical or expermintal based on the H value
used. 


%}

syms n t
alpha = 4.819e-5;

times = linspace(time_Alum(1),time_Alum(end),30); % time array
fourier_n = 20; % number  of fourier terms.

x_loc = [ 0 loc ].*0.0254 ;


for i = 1:length(times)
    
    for j = 1:length(x_loc)
    
            
% lambda_n = ((2*(k)-1)*pi)/(2*L) ;
% bn = ((-1).^(k) *8*H*L) / (( 2*(k) - 1) .^2 * pi^2);
fourier_loop_analytical = 0;
fourier_loop_exp = 0;

for f = 1:fourier_n
    
lambda_n_analytical = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_analytical = ((-1).^(f) *8*H_Alum(1)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

lambda_n_exp = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_exp = ((-1).^(f) *8*H_Alum(2)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

   fourier_loop_analytical = fourier_loop_analytical + bn_analytical.*sin( lambda_n_analytical*x_loc(j) ) * exp(- ((lambda_n_analytical)^2) *alpha * times(i)) ;
   fourier_loop_exp = fourier_loop_exp + bn_exp.*sin( lambda_n_exp*x_loc(j) ) * exp(- ((lambda_n_exp)^2) *alpha * times(i)) ;

end


     u_analytical_alum{j}{i,1} = double(T0_Alum + H_Alum(1)*x_loc(j) + fourier_loop_analytical);
     u_exp_alum{j}{i,1} = double(T0_Alum + H_Alum(2)*x_loc(j) + fourier_loop_exp);

            
    end
    
    
end
   
figure(2)

plot(times,cell2mat(u_analytical_alum{1}),'DisplayName','x0 - analytical','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_analytical_alum{n}),'DisplayName',[ 'TH' num2str(n) ' - analytical'],'LineWidth',2)

end

hold on

plot(times,cell2mat(u_exp_alum{1}),'--','DisplayName','x0 - experimental','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_exp_alum{n}),'--','DisplayName',[ 'TH' num2str(n) ' - experimental'],'LineWidth',2)

end

title('Analytical and experimental temperature variation: Aluminum')
xlabel('Time [s]')
ylabel('temperature [C^o]')


legend('Orientation','horizontal','NumColumns',3,'Location','SouthEast')
grid minor

% save plot
 set(gcf, 'Position', get(0, 'Screensize'));
 print(gcf,'figure2.png','-dpng','-r300');

 %% repeate for Brass
 
 syms n t
alpha = 35.6e-5;

times = linspace(time_Brass(1),time_Brass(end),30); % time array
fourier_n = 20; % number  of fourier terms.

x_loc = [ 0 loc ].*0.0254 ;


for i = 1:length(times)
    
    for j = 1:length(x_loc)
    
            
% lambda_n = ((2*(k)-1)*pi)/(2*L) ;
% bn = ((-1).^(k) *8*H*L) / (( 2*(k) - 1) .^2 * pi^2);
fourier_loop_analytical = 0;
fourier_loop_exp = 0;

for f = 1:fourier_n
    
lambda_n_analytical = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_analytical = ((-1).^(f) *8*H_Brass(1)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

lambda_n_exp = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_exp = ((-1).^(f) *8*H_Brass(2)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

   fourier_loop_analytical = fourier_loop_analytical + bn_analytical.*sin( lambda_n_analytical*x_loc(j) ) * exp(- ((lambda_n_analytical)^2) *alpha * times(i)) ;
   fourier_loop_exp = fourier_loop_exp + bn_exp.*sin( lambda_n_exp*x_loc(j) ) * exp(- ((lambda_n_exp)^2) *alpha * times(i)) ;

end


     u_analytical_Brass{j}{i,1} = double(T0_Brass + H_Brass(1)*x_loc(j) + fourier_loop_analytical);
     u_exp_Brass{j}{i,1} = double(T0_Brass + H_Brass(2)*x_loc(j) + fourier_loop_exp);

            
    end
    
    
end
   
figure(3)

plot(times,cell2mat(u_analytical_Brass{1}),'DisplayName','x0 - analytical','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_analytical_Brass{n}),'DisplayName',[ 'TH' num2str(n) ' - analytical'],'LineWidth',2)

end

hold on

plot(times,cell2mat(u_exp_Brass{1}),'--','DisplayName','x0 - experimental','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_exp_Brass{n}),'--','DisplayName',[ 'TH' num2str(n) ' - experimental'],'LineWidth',2)

end

title('Analytical and experimental temperature variation: Brass')
xlabel('Time [s]')
ylabel('temperature [C^o]')


legend('Orientation','horizontal','NumColumns',3,'Location','SouthEast')
grid minor
ylim([T0_Brass-4 T0_Brass+17])
% save plot
 set(gcf, 'Position', get(0, 'Screensize'));
 print(gcf,'figure3.png','-dpng','-r300');

 
 
 %% repeate for Steel
 
 syms n t
alpha = 4.05e-5;

times = linspace(time_Steel(1),time_Steel(end),30); % time array
fourier_n = 20; % number  of fourier terms.

x_loc = [ 0 loc ].*0.0254 ;


for i = 1:length(times)
    
    for j = 1:length(x_loc)
    
            
% lambda_n = ((2*(k)-1)*pi)/(2*L) ;
% bn = ((-1).^(k) *8*H*L) / (( 2*(k) - 1) .^2 * pi^2);
fourier_loop_analytical = 0;
fourier_loop_exp = 0;

for f = 1:fourier_n
    
lambda_n_analytical = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_analytical = ((-1).^(f) *8*H_Steel(1)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

lambda_n_exp = ((2*(f)-1)*pi)/(2*L.*0.0254) ;
bn_exp = ((-1).^(f) *8*H_Steel(2)*L*0.0254) / (( 2*(f) - 1) .^2 * pi^2);

   fourier_loop_analytical = fourier_loop_analytical + bn_analytical.*sin( lambda_n_analytical*x_loc(j) ) * exp(- ((lambda_n_analytical)^2) *alpha * times(i)) ;
   fourier_loop_exp = fourier_loop_exp + bn_exp.*sin( lambda_n_exp*x_loc(j) ) * exp(- ((lambda_n_exp)^2) *alpha * times(i)) ;

end


     u_analytical_Steel{j}{i,1} = double(T0_Steel + H_Steel(1)*x_loc(j) + fourier_loop_analytical);
     u_exp_Steel{j}{i,1} = double(T0_Steel + H_Steel(2)*x_loc(j) + fourier_loop_exp);

            
    end
    
    
end
   
figure(4)

plot(times,cell2mat(u_analytical_Steel{1}),'DisplayName','x0 - analytical','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_analytical_Steel{n}),'DisplayName',[ 'TH' num2str(n) ' - analytical'],'LineWidth',2)

end

hold on

plot(times,cell2mat(u_exp_Steel{1}),'--','DisplayName','x0 - experimental','LineWidth',2)
hold on
for n = 2:9
    
plot(times,cell2mat(u_exp_Steel{n}),'--','DisplayName',[ 'TH' num2str(n) ' - experimental'],'LineWidth',2)

end

title('Analytical and experimental temperature variation: Steel')
xlabel('Time [s]')
ylabel('temperature [C^o]')
ylim([12 68])


legend('Orientation','horizontal','NumColumns',3,'Location','SouthEast')
grid minor

% save plot
 set(gcf, 'Position', get(0, 'Screensize'));
 print(gcf,'figure4.png','-dpng','-r300');
