clear,clc,close all;
%% Data-analyysin harjoitustyöt - Kalteva taso
l = 6;          % toistomittausten lkm
% mittaukset

% sylinterin mittaukset

syl = 1e-3*[584.29	538.27	486.24	460.23	434.22	414.21	394.2	380.19	366.18	352.18
584.29	538.27	486.24	460.23 434.22	414.21	394.2	380.19	366.18	352.18
584.29	538.27	486.24	462.23	434.22	414.21	394.2	380.19	366.18	352.18
584.29	538.27	486.24	460.23	436.22	414.21	396.2	380.19	368.18	352.18
584.29	536.27	486.24	460.23	434.22	414.21	396.2	380.19	366.18	354.18
584.29	536.27	486.24	460.23	434.22	414.21	396.2	380.19	368.18	354.18];

% pallon mittaukset

ball = 1e-3*[590.3	550.28	490.25	462.23	438.22	416.21	396.2	382.19	368.18	354.18
586.29	536.27	494.25	460.23	436.22	416.21	396.2	382.19	368.18	352.18
584.29	534.27	492.25	460.23	436.22	414.21	396.2	378.19	370.19	356.18
584.29	534.27	490.25	466.23	438.22	418.21	398.2	380.19	368.18	354.18
582.29	534.27	488.24	460.23	438.22	414.21	400.2	382.19	368.18	354.18
584.29	534.27	490.25	462.23	440.22	418.21	394.2	382.19	366.18	354.18];

% särmiön mittaukset

sarmit = 1e-3*[686.34	594.3	512.26	500.25	458.23	448.22	408.2	390.2
668.33	594.33	534.27	496.25	458.23	432.22	406.2	398.2
686.34	588.29	536.27	498.25	456.23	438.22	412.21	394.2
682.34	594.3	536.27	502.25	460.23	438.22	412.21	402.2
700.35	590.3	542.27	504.25	454.23	438.22	416.21	404.2
686.34	600.3	542.27	496.25	458.23	436.22	414.21	400.2];

% mittausten otoskeskiarvot, otannan hajonta ja keskiarvon keskivirhe

z1 = zeros(10,1);
z2 = zeros(10,1);
v = zeros(10,2);
s = zeros(10,2);
zsar = zeros(8,1);
vsar = zeros(8,1);
ssar = zeros(8,1);

for j = 1:10
    z1(j,1) = sum(syl(:,j))/l;
    z2(j,1) = sum(ball(:,j))/l;

    s(j,1) = sqrt(sum((syl(:,j) - z1(j,1)).^2) / (l - 1));
    s(j,2) = sqrt(sum((ball(:,j) - z2(j,1)).^2) / (l - 1));
if s(j,1)==0
    s(j,1) = eps;
end
if s(j,2)==0
    s(j,2) = eps;
end
    v(j,1) = s(j,1)/(sqrt(l));
    v(j,2) = s(j,2)/(sqrt(l));

end
for j = 1:8
    zsar(j) = sum(sarmit(:,j))/l;

    ssar(j) = sqrt(sum((sarmit(:,j) - zsar(j)).^2) / (l - 1));

    vsar(j) = ssar(j)/(sqrt(l));
end

% painomatriisi W yleistetylle Gauss-Newtonille
W = zeros(10);
WW = zeros(10);
Wsar = zeros(8);
for j = 1:10
    for ii = 1:10
        W(j,j) = 1/s(j,1)^2;
        WW(j,j) = 1/s(j,2)^2;
    end
end
for j = 1:8
    for ii = 1:8
        Wsar(j,j) = 1/ssar(j)^2;
    end
end
sylball = [10	584.29	585.2916667
12	537.6033333	537.2716667
14	486.24	490.915
16	460.5633333	461.8966667
18	434.5533333	437.8866667
20	414.21	416.21
22	395.2	396.8666667
24	380.19	381.19
26	366.8466667	368.1816667
28	352.8466667	354.18
];

sylball(:,2) = 1e-3*sylball(:,2);
sylball(:,3) = 1e-3*sylball(:,3);

% särmiön mittaukset
sar = [20	685.066667
22	593.6366667
24	533.935
26	499.5833333
28	457.5633333
30	438.5533333
32	411.54
34	398.2
];

sar(:,2) = 1e-3*sar(:,2);

% mittausten matriisit taulukoksi käsittelyä varten

%meas1 = array2table(sylball,'VariableNames',{'kulma','syl','ball'});
%meas2 = array2table(sar,'VariableNames',{'kulma','sar'});

% sylinterin dimensiot

% h = 7e-2;          % sylinterin korkeus (m)
dsyl = 4.76e-2;    % sylinterin halkaisija (m)
rsyl = dsyl/2;  % sylinterin säde (m)
msyl = 119e-3;     % sylinterin massa (kg)


% pallon dimensiot

dball = 6.3e-2;    % pallon halkaisija (m)
rball = dball/2;    % pallon säde (m)
mball = 190e-3;    % pallon massa (kg)


% suorakulmion dimensiot

hh = 3e-2;         % suorakulmion korkeus (m)
ll = 9.55e-2;      % suorakulmion pituus (m)
ww = 7.38e-2;      % suorakulmion leveys (m)
mm = 199e-3;       % suorakulmion massa (kg)

% teoreettiset hitausmomentit

Ith = 1/2*msyl*rsyl^2;          % sylinterin teoreettinen hitausmomentti
IIth = 2/5*mball*rball^2;       % pallon teoreettinen hitausmomentti

% teoreettisten hitausmomenttien virhe
vm = 0.0005;        % vaa'asta johtuva virhe
vr = 0.00005;       % työntömitasta johtuva virhe

% sylinterin teor hitausmomentin virhe
vith = sqrt((1/2*rsyl^2)^2*vm^2 + (msyl*rsyl)^2*vr^2);
% pallon teor hitausmomentin virhe
viith = sqrt((2/5*rball^2)^2*vm^2 + (4/5*mball*rball)^2*vr^2);

%----------------------------------------


% vakiot

g = 9.81;           % putoamiskiihtyvyys maan pinnalla

s1 = 0.13;         % valoportin 1 etäisyys lähtöpisteestä (m)
s2 = 0.575;         % valoportin 2 etäisyys lähtöpisteestä (m)

Csyl = sqrt(2) * (sqrt(s2 - dsyl) - sqrt(s1 - dsyl));
Cball = sqrt(2) * (sqrt(s2 - dball) - sqrt(s1 - dball));
Csar = sqrt(2) * (sqrt(s2 - ll) - sqrt(s1 - ll));

%====================================================

n = 500;             % iterointien enimmäismäärä
%ask = 0.03;            % askelparametri

figure(1)
tiledlayout(2,1)

%% Gauss-Newtonin iterointi sylinterille

% Alkuarvaus hitausmomentille
I = 0.00015;
Isyl = I;
max_iter = 500;
ask = 0.3;


iterlog = [];
Iloki = [];

for j = 1:max_iter
    % Havaintomalli
    alpha = sind(sylball(:,1));
    h1 = Csyl * sqrt(1 + Isyl / (msyl * rsyl^2)) ./ sqrt(g * alpha);
    
    % Jacobin matriisi
    J = Csyl ./ (2 * sqrt(g * alpha)) .* (1 / (msyl * rsyl^2)) ./ sqrt(1 + Isyl / (msyl * rsyl^2));
    
    % Residuaali
    r = syl(1,:).' - h1;
    
    % Gauss-Newton päivitys
    deltaI = inv(J.' * J) * J.' * r;
    Isyl = Isyl + ask * deltaI;
    
    % tallennetaan iteraatiot
    iterlog(end+1) = j;
    Iloki(end+1) = Isyl;
    
    % suppenemisehto
    if abs(deltaI) < 1e-8
        break
    end
end


nexttile
plot(iterlog, Iloki, '-o', 'LineWidth', 1.5), hold on, grid on
yline(Ith, 'LineWidth', 3)
title('Gauss-Newton sylinterille','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraation arvo (I_{cm})','FontSize',14)
legend('GN menetelmän iteraatiot','Teoreettinen hitausmomentti',location='best')

% Kovarianssi ja virhearvio
Vsyl = (norm(z1-h1))/(10-1);

kov11 = Vsyl*inv(J.'*J);
sylvirhe = sqrt(kov11);


%% Gauss-Newtonin iterointi sylinterille (TÄSSÄ ON VIRHE, ÄLÄ KÄYTÄ)

% figure(1)
% tiledlayout(2,1)
% % käytetään jokaisen mittaussarjan ensimmäistä mittausta
% 
% % alustetaan sylinterille iteraatiovektori
% Isyl = [];
% Isyl(1) = I;
% j1 = [];
% J1save = [];
% 
% for j=1:n
%     h1 = Csyl * (1 + Isyl(j) / (msyl * rsyl^2)) ./ (sqrt(g * sind(sylball(:,1))));
% 
%     J1 = Csyl ./ (2 * sqrt(g * sind(sylball(:,1)))) .* (1 / (msyl * rsyl^2)) ./ sqrt(1 + Isyl(j)/(msyl * rsyl^2));
% 
%     Isyl(j+1) = Isyl(j) + ask * inv(J1.' * J1) * J1.' * (syl(1,:).' - h1);
% 
%     j1 = [j1;j];
% 
%     J1save = [J1save,J1];
% 
%     if round(Isyl(j+1),7)==round(Isyl(j),7)
%         break
%     end
% end
% j1 = [j1;j+1];
% 
% nexttile
% plot(j1,Isyl,'-o','LineWidth',1.5)
% hold on, grid on
% plot(j1,Ith*ones(j+1,1),'LineWidth',3)
% title('Gauss-Newton sylinterille','FontSize',16,'FontWeight','bold')
% xlabel('Iteraatioiden lkm','FontSize',14)
% ylabel('Iteraation arvo (I_{cm})','FontSize',14)
% legend('GN menetelmän iteraatiot','Teoreettinen hitausmomentti',location='best')
% 
% 
% % estimaatin virhe (kovarianssi)
% % Rsyl = diag(sylball(:,2).^2);
% 
% Vsyl = (norm(z1-h1))/(10-1);
% 
% kov11 = Vsyl*inv(J1save(:,20).'*J1save(:,20));
% sylvirhe = sqrt(kov11);
%% Yleistetty Gauss-Newton sylinterille

IIlog = [];
IIlog(1) = I;
IIsyl = I;
jj1 = [];
JJ1save = [];
for jj = 1:n
    % Mallifunktio h(I)
    alpha = sind(sylball(:,1));
    hh1 = Csyl * sqrt(1 + IIsyl / (msyl * rsyl^2)) ./ sqrt(g * alpha);
    
    % Jacobin matriisi dh/dI (10x1)
    JJ1 = Csyl ./ (2 * sqrt(g * alpha)) .* (1 / (msyl * rsyl^2)) ./ sqrt(1 + IIsyl / (msyl * rsyl^2));
    
    % Residuaali
    r = z1 - hh1;
    
    % Gauss-Newton päivitys
    deltaII = inv(JJ1.' * JJ1) * JJ1.' * r;
    IIsyl = IIsyl + ask * deltaII;
    
    % Lokitietoa
    jj1(end+1) = jj;
    IIlog(end+1) = IIsyl;
    
    % suppenemisehto
    if abs(deltaII) < 1e-8
        break
    end
end
jj1(end+1) = jj+1;

nexttile
plot(jj1,IIlog,'-o','LineWidth',1.5)
hold on, grid on
yline(Ith,'LineWidth',3)
title('Yleistetty Gauss-Newton sylinterille','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraation arvo (I_{cm})','FontSize',14)
legend('YGN menetelmän iteraatiot','Teoreettinen hitausmomentti',location='best')


% estimaatin virhe (kovarianssi)

kov12 = inv(JJ1.' * W * JJ1);

sylvirheyl = sqrt(kov12);
%% Linearisointi sylinterille (PNS, yleistetty PNS)
figure(2)
tiledlayout(2,1)

z = (syl(1,:).^(-1)).^2.';

x = sind(sylball(:,1));
A = inv(x.' * x) *x.' * z;

AA = inv(x.' * W * x)* x.'* W * (z1.^(-1)).^2;

Jsyl = (g/(Csyl^2*A) - 1)*(msyl*rsyl^2);
JJsyl = (g/(Csyl^2*AA) - 1)*(msyl*rsyl^2);


ismo = fitlm(A*x,z);
seppo = fitlm(AA*x,z);

nexttile
plot(ismo)
grid on
ylabel('1/t^2','FontSize',14)
xlabel('sin(alpha)','FontSize',14)
title('Linearisoidulla PNS-menetelmällä','FontSize',16,'FontWeight','bold')



nexttile
plot(seppo)
grid on
xlabel('sin(alpha)','FontSize',14)
ylabel('1/t^2','FontSize',14)
title('Linearisoidulla yleistetyllä PNS-menetelmällä','FontSize',16,'FontWeight','bold')


% estimaatin virhe (kovarianssi)
Vsyl1 = (norm(z-A*x))/(10-1);
dJdA = -g*msyl*rsyl^2/(Csyl^2*A^2);
dJdAA = -g*msyl*rsyl^2/(Csyl^2*AA^2);

kov13 = Vsyl1*inv(x.' * x);
kov14 = inv(x.' * W * x);

syllinvirhe = sqrt(kov13);
J0virhelin = sqrt(dJdA^2*syllinvirhe^2);

syllinvirheyl = sqrt(kov14);
J0virhelinyl = sqrt(dJdAA^2*syllinvirheyl^2);
%% Gauss-Newton pallolle
figure(3)
tiledlayout(2,1)
i1 = [];

Iball = 0.000015;
Iballok = Iball;


ask = 0.03;
for i = 1:n
    % Havaintofunktio h
    alpha = sind(sylball(:,1));
    h2 = Cball * sqrt(1 + Iball / (mball * rball^2)) ./ sqrt(g * alpha);

    % Jacobin matriisi
    J2 = Cball ./ (2 * sqrt(g * alpha)) .* (1 / (mball * rball^2)) ./ sqrt(1 + Iball / (mball * rball^2)).';

    % residuaali
    resball = ball(1,:).' - h2;

    % Gauss-Newtonin iteraatio
    deltaIball = inv(J2.' * J2) * J2.' * resball;
    Iball = Iball + ask * deltaIball;

    Iballok(end+1) = Iball;
    i1 = [i1;i];

    if abs(deltaIball) < 1e-8
        break
    end
end
i1 = [i1;i+1];

nexttile
plot(i1,Iballok,'-o','LineWidth',1.5)
hold on,grid on
yline(IIth,'LineWidth',3)
title('Gauss-Newton pallolle','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraation arvo (I_{cm})','FontSize',14)
legend('GN menetelmän iteraatiot','Teoreettinen hitausmomentti',location='best')


% estimaatin virhe (kovarianssi)
Rball = diag(sylball(:,3).^2);

Vball = (norm(z2 - h2))/(10-1);

kov21 = Vball * inv(J2.'*J2);

ballvirhe = sqrt(kov21);
%% Yleistetty Gauss-Newton pallolle
ii1 = [];
IIball = 0.000015;
IIballok = IIball;


JJ2save = [];

for ii = 1:n
    % Havaintomalli
    h2 = Cball * sqrt(1 + IIball / (mball * rball^2)) ./ (sqrt(g * sind(sylball(:,1))));

    % Jacobin matriisi
    JJ2 = Cball ./ (2 * sqrt(g * sind(sylball(:,1)))) .* 1 / (mball * rball^2) * 1/(sqrt(1 + IIball / (mball * rball^2))).';

    % Residuaali
    resball2 = z2 - h2;

    % Yleistetty Gauss-Newtonin iteraatio
    deltaIIball = inv(JJ2.' * WW * JJ2) * JJ2.' * WW * resball2;
    IIball = IIball + ask * deltaIIball;

    IIballok(end + 1) = IIball;
    ii1 = [ii1;ii];

    JJ2save = [JJ2save,JJ2];
    if abs(deltaIIball) < 1e-8
        break
    end
end

ii1 = [ii1;ii+1];

nexttile
plot(ii1,IIballok,'-o','LineWidth',1.5)
hold on,grid on
yline(IIth,'LineWidth',3)
title('Yleistetty Gauss-Newton pallolle','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraation arvo (I_{cm})','FontSize',14)
legend('YGN menetelmän iteraatiot','Teoreettinen hitausmomentti',location='best')


% estimaatin virhe (kovarianssi)

kov22 = inv(JJ2.' * WW * JJ2);

ballvirheyl = sqrt(kov22);

%% Linearisointi pallolle (PNS, yleistetty PNS)
figure(4)
tiledlayout(2,1)

zz = (ball(1,:).^(-1)).^2.';

% A = g/(Cball^2*(1+J0/(mr^2)));
x = sind(sylball(:,1));
A = inv(x.' * x)*x.' * zz;

AA = inv(x.' * WW * x)*x.' * WW * (z2.^(-1)).^2;

Jball = (g/(Cball^2*A)-1)*mball*rball^2;
JJball = (g/(Cball^2*AA)-1)*mball*rball^2;

ossi = fitlm(A*x,zz);
kalle = fitlm(AA*x,zz);

nexttile
plot(ossi)
grid on
xlabel('sin(alpha)','FontSize',14)
ylabel('1/t^2','FontSize',14)
title('Linearisoitu PNS-menetelmä pallolle','FontSize',16,'FontWeight','bold')

nexttile
plot(kalle)
grid on
xlabel('sin(alpha)','FontSize',14)
ylabel('1/t^2','FontSize',14)
title('Linearisoitu yleistetty PNS-menetelmä pallolle','FontSize',16,'FontWeight','bold')


% estimaatin virhe (kovarianssi)
Vball1 = (norm(zz - A*x))/(10-1);
dJdA = -g*mball*rball^2/(Cball^2*A^2);
dJdAA = -g*mball*rball^2/(Cball^2*AA^2);

kov23 = Vball1*inv(x.' * x);
kov24 = inv(x.' * WW * x);

ballinvirhe = sqrt(kov23);
J0ballvirhe = sqrt(dJdA^2*ballinvirhe^2);

ballinvirheyl = sqrt(kov24);
J0ballvirheyl = sqrt(dJdAA^2*ballinvirhe^2);
%% Särmiön Gauss-Newtonin iteraatiot
figure(5)
tiledlayout(2,1)
k1 = [];
mu0 = 0.1;           % alkuarvaus särmiön liikekitkakertoimelle
mu = [mu0];             % alustetaan vektori liikekitkakertoimen iteraatioille

J3save = [];
ask1 = 0.1;
for k = 1:n
    h3 = Csar / sqrt(g) * 1 ./ sqrt(sind(sar(:,1)) - mu(k) * cosd(sar(:,1)));
    J3 = Csar / (2*sqrt(g)) .* cosd(sar(:,1)) ./ (sind(sar(:,1)) - mu(k) .* cosd(sar(:,1)) .* sqrt(sind(sar(:,1)) - mu(k) .* cosd(sar(:,1))));
    mu(k+1) = mu(k) + ask1 * inv(J3.' * J3) * J3.' * (sarmit(1,:).' - h3);

    k1 = [k1;k];

    J3save = [J3save,J3];
    if round(mu(k+1),4) == round(mu(k),4)
        break
    end
end

k1 = [k1;k+1];

nexttile
plot(k1,mu,'-o','LineWidth',1.5)
grid on
title('Gauss-Newton särmiölle','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraatioiden arvo (\mu)','FontSize',14)


% estimaatin virhe (kovarianssi)
Rsar = diag(sar(:,2).^2);
Vsar = (norm(zsar - h3))^2/(8-1);

kov31 = Vsar*inv(J3.' * J3);


sarvirhe = sqrt(kov31);
%% Yleistetty Gauss-Newton särmiölle

kk1 = [];
muu = [mu0];

JJ3save = [];

for kk = 1:n
    h3 = Csar / sqrt(g) * 1 ./ sqrt(sind(sar(:,1)) - muu(kk) * cosd(sar(:,1)));
    JJ3 = Csar / (2*sqrt(g)) .* cosd(sar(:,1)) ./ (sind(sar(:,1)) - muu(kk) .* cosd(sar(:,1)) .* sqrt(sind(sar(:,1)) - muu(kk) .* cosd(sar(:,1))));
    muu(kk+1) = muu(kk) + ask1 * inv(JJ3.'* Wsar * JJ3) * JJ3.'* Wsar * (zsar - h3);

    kk1 = [kk1;kk];

    JJ3save = [JJ3save,JJ3];
    if round(muu(kk+1),4) == round(muu(kk),4)
        break
    end
end

kk1 = [kk1;kk+1];


nexttile
plot(kk1,muu,'-o','LineWidth',1.5)
grid on
title('Yleistetty Gauss-Newton särmiölle','FontSize',16,'FontWeight','bold')
xlabel('Iteraatioiden lkm','FontSize',14)
ylabel('Iteraatioiden arvo (\mu)','FontSize',14)

% estimaatin virhe (kovarianssi)

kov32 = inv(JJ3save(:,15).' * Wsar * JJ3save(:,15));

sarvirheyl = sqrt(kov32);