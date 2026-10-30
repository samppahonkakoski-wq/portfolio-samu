%% ETH Project work

%% Task 1: Locating an earthquake
clc; clear; close all

% seismic measurement stations measure seismic waves at the times

t = [3.12 3.26 2.98 3.12 2.84 2.98]';   % time in s

% Covariance of noise
Cv = 0.1^2 * eye(6);

% locations of each measurement stations
loc_xy = [3 15; % meas station #1
         3 16;  % meas station #2
         4 15;  % meas station #3
         4 16;  % meas station #4
         5 15;  % meas station #5
         5 16]; % meas station #6

% speed of seismic wave
v = 5; % km/s

%% (a) Derive non-linear observation model f the problem for estimating of
% the location of the earthquake
%% (b) Determine ML estimate for the location of the earthquake

% initial guess of the location
x0 = 0;
y0 = 0;
th = [x0 y0]';





% Gauss-Newton iteration
%number of iterations
N = 7;

% initialize estimation vector
ths = zeros(2,N);
ths(:,1) = th;


% Gauss-Newton iteration
for j = 2:N     % number of measurements
    % Model function
    dx = th(1)*ones(6,1) - loc_xy(:,1);
    dy = th(2)*ones(6,1) - loc_xy(:,2);
    h = sqrt(dx.^2 + dy.^2) / v;

    % Jacobian matrix
    J = 1/v * [dx ./ sqrt(dx.^2 + dy.^2), dy ./ sqrt(dx.^2 + dy.^2)];
    
    % Residual
    r = t - h;

    % Gauss-Newton iteration
    th = th + inv(J' * inv(Cv)* J)* J' * inv(Cv) * r;

    % save the new theta value into ths vector

    ths(:,j) = th;
end

% final form of iterated model function
h = sqrt(dx.^2 + dy.^2) / v;

%% (c) Plot likelihood density with given observations


xg = linspace(min(loc_xy(:,1))-5, th(1)+6, 1000);
yg = linspace(max(loc_xy(:,2))+5, th(2)-6, 1000);
[X, Y] = meshgrid(xg, yg);

L = zeros(size(X));  % likelihood grid

for k=1:numel(L)
    thc = [X(k); Y(k)];
    dx = thc(1) - loc_xy(:,1);
    dy = thc(2) - loc_xy(:,2);
    hc = sqrt(dx.^2 + dy.^2)/v;
    L(k) = exp(-1/2*(t-hc)'*inv(Cv)*(t-hc));
end



figure(1)
contour(xg,yg,L)
hold on; grid on
plot(loc_xy(:,1), loc_xy(:,2), 'x','LineWidth',3)
plot(th(1),th(2),'o','LineWidth',3)
xlabel('x (km)')
ylabel('y (km)')
legend('Uskottavuusfunktio','Seismometrien paikat','ML-estimaatti',location='best')

%% (d) Likelihood density with one observation
% Let's choose the first observation and simulate the results with it


xg = linspace(loc_xy(1,1)+18, loc_xy(1,1)-18 , 100);
yg = linspace(loc_xy(1,2)-18 , loc_xy(1,2)+18, 100);
[X, Y] = meshgrid(xg, yg);

L = zeros(size(X));  % likelihood grid

for k=1:numel(L)
    thc = [X(k); Y(k)];
    dx = thc(1) - loc_xy(1,1);
    dy = thc(2) - loc_xy(1,2);
    hc = sqrt(dx.^2 + dy.^2)/v;
    L(k) = exp(-1/2*(t(1)-hc)'* 0.1^(-2) *(t(1)-hc));
end

figure(2)
contour(xg,yg,L)
hold on; grid on
plot(loc_xy(1,1), loc_xy(1,2),'x', 'LineWidth',3)
plot(th(1), th(2), 'o', 'LineWidth',3)
xlabel('x (km)'); ylabel('y (km)')
legend('Uskottavuusfunktio','Seismometrin paikka','ML-estimaatti 6 seismometrillä', location='best')


%% (e) MAP estimate of earthquake with prior information

% Prior mean
eta0 = [16 6]';
Cth = [1 0.7; 0.7 1];

th = eta0;

% First compute MAP estimate of earthquake location
tht = zeros(N+1,2);
tht(1,:) = th;

for i=1:N
    dx = tht(i,1) - loc_xy(:,1);
    dy = tht(i,2) - loc_xy(:,2);
    h = sqrt(dx.^2 + dy.^2)/v;

    % Jacobian
    J = 1/v * [dx ./ sqrt(dx.^2 + dy.^2), dy ./ sqrt(dx.^2 + dy.^2)];

    % residual
    r = t - h;
    th = th + inv(J' * inv(Cv) * J + inv(Cth)) * (J' * inv(Cv) * r - inv(Cth)*(th - eta0));
    tht(i+1,:) = th;

end

xg = linspace(min(loc_xy(:,1))-5, th(1)+5, 100);
yg = linspace(max(loc_xy(:,2))+5, th(2)-5, 100);
[X, Y] = meshgrid(xg, yg);

L = zeros(size(X));  % likelihood grid

for k=1:numel(L)
    thc = [X(k); Y(k)];
    dx = thc(1) - loc_xy(:,1);
    dy = thc(2) - loc_xy(:,2);
    hc = sqrt(dx.^2 + dy.^2)/v;
    L(k) = exp(-1/2*(t-h)'*inv(Cv)*(t-h)) * exp(-1/2 * (thc-eta0)' * inv(Cth) * (thc-eta0));
end

figure(3)
contour(xg,yg,L)
hold on; grid on
plot(loc_xy(:,1), loc_xy(:,2), 'x', 'LineWidth',3)
plot(th(1),th(2),'o','LineWidth',3)
xlabel('x (km)')
ylabel('y (km)')
legend('Posterioritiheys','Seismometrien paikat','MAP-estimaatti',location='best')


%% (f) Try different number of measurements and prior covariance and prior mean

% Prior mean
eta0 = [14 9]';
Cth = [4 2.2; 2.2 5];

th = eta0;

% First compute MAP estimate of earthquake location
N = 15;
tht = zeros(N+1,2);
tht(1,:) = th;

for i=1:N
    dx = tht(i,1) - loc_xy(1:5,1);
    dy = tht(i,2) - loc_xy(1:5,2);
    h = sqrt(dx.^2 + dy.^2)/v;

    % Jacobian
    J = 1/v * [dx ./ sqrt(dx.^2 + dy.^2), dy ./ sqrt(dx.^2 + dy.^2)];

    % residual
    r = t(1:5) - h;

    % Gauss-Newton
    th = th + inv(J' * inv(Cv(1:5,1:5)) * J + inv(Cth)) * (J' * inv(Cv(1:5,1:5)) * r - inv(Cth)*(th - eta0));
    tht(i+1,:) = th;

end

xg = linspace(min(loc_xy(:,1))-5, th(1)+5, 100);
yg = linspace(max(loc_xy(:,2))+5, th(2)-5, 100);
[X, Y] = meshgrid(xg, yg);

L = zeros(size(X));  % posterior grid

for k=1:numel(L)
    thc = [X(k); Y(k)];
    dx = thc(1) - loc_xy(1:5,1);
    dy = thc(2) - loc_xy(1:5,2);
    hc = sqrt(dx.^2 + dy.^2)/v;
    L(k) = exp(-1/2 * (thc-eta0)' * inv(Cth) * (thc-eta0)) * exp(-1/2 * ((t(1:5)-h)' * inv(Cv(1:5,1:5)) * (t(1:5)-h)));
end

figure(4)
contour(xg,yg,L)
hold on; grid on
plot(loc_xy(1:5,1), loc_xy(1:5,2), 'x', 'LineWidth',3)
plot(th(1),th(2),'o','LineWidth',3)
xlabel('x (km)')
ylabel('y (km)')
legend('Posterioritiheys','Seismometrien paikat','MAP-estimaatti',location='best')

