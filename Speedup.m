%% This code is intended compare "Consecutive" and "Adaptive" schemes
% Analytical solution (1D C-D equation without adsorption) is taken from Lake et al. (2014)
% Porosity is implict in pore volume - ALphi
% Units in [], example in ()
clear
%% INPUT
Ci =  1000;                                             % Initial concentration [ppm] (1000)
Cinj = 500;                                             % Injection concentration [ppm] (500)
tstep = 0.005;                                          % Timestep [PV] (0.005)                    
q = 1; ALphi=200;                                       % Injection rate [m3/day] (1) and pore volume [m3] (200)
Npe = 100;                                              % Peclet number (50-1000)
tol = 1;                                                % Tolerance for criteria evaluation [%] (5)
%% Indexing
t = linspace(0,(1/tstep)-1,1/tstep);                    % Time [day] (e.g. 1PV/0.005=200)  
tola = ones(100,size(t,2))*tol;                         % Indexing tolerance for each gridblock
%% Injected PV calculation
td = zeros(1,size(t,2));
for i=2:size(t,2)
    td(i) = q*(t(i)-t(i-1))/(ALphi) + td(i-1);
end
%% Analytical solution for concentration
xd = linspace(0,1,100);                                 % Dimensionaless distance
for i=1:size(xd,2)
    for j=1:size(t,2)
        Cd(i,j) = 0.5*erfc((xd(i)-td(j))/(2*sqrt(td(j)/Npe))) + (exp(xd(i)*Npe)/2)*erfc((xd(i)+td(j))/(2*sqrt(td(j)/Npe)));
    end
end
C = Cd*(Cinj-Ci)+Ci;
%% Normalized concentration slope calculation (omega)
for j=1:size(xd,2)
    for i=2:size(t,2)
        dC(j,i-1) = C(j,i)-C(j,i-1);
        dt(j,i-1) = t(i)-t(i-1);
        slope(j,i-1) = 100*(dC(j,i-1)/dt(j,i-1));
        norm_slope(j,i-1) = slope(j,i-1)/C(j,i);
    end
end
%% Correction factor calculation (fac)
for j=1:size(xd,2)
    for i=2:size(t,2)-1
        fac(j,i-1) = (1+abs(norm_slope(j,i)))/(1+abs(norm_slope(j,i-1)));
        tola(j,i) = tola(j,i-1)/(fac(j,i-1));
    end
end
%% Speedup scheme calculation and criteria evaluation
Cs1 = ones(size(xd,2),1)*Ci;
Cs2 = ones(size(xd,2),1)*Ci;
Call1 = nan(size(xd,2),size(t,2));
Call2 = nan(size(xd,2),size(t,2));
Call3 = nan(size(xd,2),size(t,2));
for j=1:size(xd,2)
    for i=2:size(t,2)-1
        error1(j,i)=100*abs(((C(j,i)-C(j,i-1))/C(j,i-1)));  % Consecutive
        if error1(j,i)>tol
            Call1(j,i)=C(j,i);
        end
        error2(j,i)=100*abs(((C(j,i)-Cs1(j))/Cs1(j)));      % Cumulative
        if error2(j,i)>tol
            Call2(j,i)=C(j,i);
            Cs1(j)=C(j,i);
        end
        error3(j,i)=100*abs(((C(j,i)-Cs2(j))/Cs2(j)));      % Adaptive
        if error3(j,i)>tola(j,i)         
            Call3(j,i)=C(j,i);
            Cs2(j)=C(j,i);
        end
    end
end
%% Number of calls for each scheme
num_calls1 = sum(~isnan(Call1(50,:)));                      % Consecutive
num_calls2 = sum(~isnan(Call2(50,:)));                      % Cumulative
num_calls3 = sum(~isnan(Call3(50,:)));                      % Adaptive
%% Plots
% Consecutive
figure (1);
plot(td,C(size(xd,2)/2,:),'LineWidth',2)
hold on
ylabel('Concentration [ppm]')
xlabel('Pore volume [PV]')
plot(td,Call1(size(xd,2)/2,:),'k.','MarkerSize',15)
legend('C-D solution','error>tol','Location','southwest')
title('Concentration profile at xd=0.5 - Consecutive')
% Cumulative
figure (2);
plot(td,C(size(xd,2)/2,:),'LineWidth',2)
ylabel('Concentration [ppm]')
xlabel('Pore volume [PV]')
hold on
plot(td,Call2(size(xd,2)/2,:),'k.','MarkerSize',15)
legend('C-D solution','error>tol','Location','southwest')
title('Concentration profile at xd=0.5 - Cumulative')
% Adaptive
figure (3);
yyaxis left
plot(td,C(size(xd,2)/2,:),'LineWidth',2)
hold on
plot(td,Call3(size(xd,2)/2,:),'k.','MarkerSize',15)
ylabel('Concentration [ppm]')
xlabel('Pore volume [PV]')
yyaxis right
plot(td,tola(size(xd,2)/2,:),'LineWidth',2)
legend('C-D solution','error>tol','tol','Location','southwest')
title('Concentration profile at xd=0.5 - Adaptive')
ylabel('Tolerance [%]')
% Number of calls comparison
figure (4)
X = categorical({'Consecutive','Cumulative','Adaptive'});
X = reordercats(X,{'Consecutive','Cumulative','Adaptive'});
Y = [num_calls1 num_calls2 num_calls3];
bar(X,Y,0.4)
ylabel('Number of calls')

