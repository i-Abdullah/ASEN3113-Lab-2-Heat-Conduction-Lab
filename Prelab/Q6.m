%% prelab 2: ASEN 3113 Thermodynamics


%% hoeskeeping

clear
clc
close all

%% problem 3:

th1_x = 1+(3/8); % location of first tehrmocouple:


T = [ 18.53 22.47 26.87 30.05 35.87 38.56 41.50 46.26 ];
x = [ th1_x:0.5:th1_x+(0.5*7) ];


% linear fit:


[xData, yData] = prepareCurveData( x, T );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Temp vs. tehrmocouples', 'Linear fit', 'Location', 'NorthWest', 'Interpreter', 'none' );
title('thermocouple measurements for steady state')
% Label axes
xlabel( 'location of thermocouple [in]', 'Interpreter', 'none' );
ylabel( 'Temp [C]', 'Interpreter', 'none' );
grid on

T0 = fitresult.p2;
H = fitresult.p1;
H = H/0.0254;

%% question 5

L = x(end) + 1 ; % inches
x = L - 1; 

L = L*0.0254 ;
x = x * 0.0254;

syms n t
bn = ((-1)^n *8*H*L) / (( 2*n - 1) ^2 * pi^2);
lambda_n = ((2*n-1)*pi)/2*L ;
%alpha = 8e-5;

%f = bn*sin(lambda_n*x)*exp(-(lambda_n)^2*alpha*t)

%f = (((-1)^n *8*H*L) / (( 2*n - 1) ^2 * pi^2)) * sin( ((2*n-1)*pi)/2*L*x)*exp(- ( ((2*n-1)*pi)/2*L)^2*alpha * t) ;

n_array = [ 1 ];
t = [1:2:1000];


alpha_values = [ 1:5:90 ] .* 10^-5 ; 
%alpha_values = [8e-5];
u = zeros(length(n_array),length(t));

for z = 1:length(alpha_values);
    
    alpha = alpha_values(z);
    
    
    
for i = 1:length(t)
    
    
    for j = 1:length(n_array)
        
        sum_2 = 0 ;
        
        for k = 1:n_array(j)
            
            
lambda_n = ((2*(k)-1)*pi)/(2*L) ;
bn = ((-1).^(k) *8*H*L) / (( 2*(k) - 1) .^2 * pi^2);





sum_2 = sum_2 +  bn.*sin( lambda_n*x ) *exp(- ( lambda_n )^2*alpha * t(i)) ;
            
        end
        
        
    end
  
       u(j,i) = T0 + H*x + sum_2;

            
end


plot(t,u,'DisplayName',['alpha = ' num2str(alpha)],'LineWidth',2);
hold on


end

legend
grid minor
xlabel('time [s]')
ylabel('T [C]')
title('Temprature profile with time at the location of the last thermocouple for thermal diffusivity')


%% plots



% figure(2)
% plot(alpha_values(:),u(:,1),'-*r','LineWidth',2)
% xlabel('\alpha')
% ylabel('T [C]')
% title('t = 1s')
% grid minor
% hold on
% 
% figure(3)
% plot(alpha_values,u(:,2),'-*b','LineWidth',2)
% xlabel('\alpha')
% ylabel('T [C]')
% title('t = 1000s')
% grid minor
% 
% 
% hold on
% 
