clear all
close all

mfreq=1550;
mtheta=31;
fmin= 0.04;
fmax= 0.25;
h_gen=9;
delta=2;
fm = 0.1;
theta_input=00;
sigma_theta_input=1;
gamma_spec=3;
DY = 1;
Nglob=2999;
alpha_c=0;
PERIODIC=1;
SMALL=eps('single');

[Freq0,theta0,Hmo_each0]=generate_salitin_freq_dire(mfreq,mtheta,fmin,fmax,h_gen,delta,fm,theta_input,sigma_theta_input,gamma_spec,DY,Nglob,alpha_c,PERIODIC,SMALL);


[Freq,theta,Hmo_each,Freq_new,theta_new,Hmo_each_new]=generate_modified_salitin_freq_dire(mfreq,mtheta,fmin,fmax,h_gen,delta,fm,theta_input,sigma_theta_input,gamma_spec,DY,Nglob,alpha_c,PERIODIC,SMALL);


% interpolate each to estimate directional shape at peak frequency
td=[-60:1:60];
t =td*pi/180;
f = fm+0*t;

F0 = scatteredInterpolant(Freq0',theta0',Hmo_each0'.^2);
E0 = F0(f,t);

F1 = scatteredInterpolant(Freq',theta',Hmo_each'.^2);
E1 = F1(f,t);

F2 = scatteredInterpolant(Freq_new',theta_new',Hmo_each_new'.^2);
E2 = F2(f,t);

fig0 = figure;
p0 = plot(td,E0,'-k',td,E2,'-r','linewidth',2);

% fit cos^2s(theta)
ft = fittype( 'A*(cos(x + C)^2)^s', ...
              'independent','x', ...
              'coefficients',{'A','C','s'});

fit0 = fit(t',E0', ft,'StartPoint', [max(E0)-min(E0), 0, 1]);
fit2 = fit(t',E2', ft,'StartPoint', [max(E2)-min(E2), 0, 1]);

hold on,
p1 = plot(td,fit0(t),'--k',td,fit2(t),'--r');


spreads1 = [sqrt(2/(fit0.s+1))*(180/pi), sqrt(2/(fit2.s+1))*(180/pi)]


df = mtheta*(fmax-fmin)/(mfreq-1);

i0 = find(abs(Freq0-fm)<2*df);
i1 = find(abs(Freq-fm)<2*df);
i2 = find(abs(Freq_new-fm)<2*df);

% $$$ m0 = sum( theta0(i0).*Hmo_each0(i0).^2 )/numel(i0);
% $$$ s0 = sqrt(sum( (theta0(i0)-m0).^2.*Hmo_each0(i0).^2 )/(numel(i0)-1));
% $$$ 
% $$$ m1 = sum( theta(i1).*Hmo_each(i1).^2 )/numel(i1);
% $$$ s1 = sqrt(sum( (theta(i1)-m1).^2.*Hmo_each(i1).^2 )/(numel(i1)-1));
% $$$ 
% $$$ m2 = sum( theta_new(i2).*Hmo_each_new(i2).^2 )/numel(i2);
% $$$ s2 = sqrt(sum( (theta_new(i2)-m2).^2.*Hmo_each_new(i2).^2 )/(numel(i2)-1));
% $$$ 
% $$$ [s0 s1 s2]*180/pi



fitresult0 = fit(theta0(i0)',Hmo_each0(i0)'.^2, ft,'StartPoint', [max(Hmo_each0(i0).^2)-min(Hmo_each0(i0).^2), 0, 1]);
fitresult1 = fit(theta(i1)',Hmo_each(i1)'.^2, ft,'StartPoint', [max(Hmo_each(i1).^2)-min(Hmo_each(i1).^2), 0, 1]);
fitresult2 = fit(theta_new(i2)',Hmo_each_new(i2)'.^2, ft,'StartPoint', [max(Hmo_each_new(i2).^2)-min(Hmo_each_new(i2).^2), 0, 1]);


figure, plot(theta0(i0)*180/pi,Hmo_each0(i0).^2,'*k',theta_new(i2)*180/pi,Hmo_each_new(i2).^2,'or')
hold on,plot(td,fitresult0(t),'-k',td,fitresult2(t),'-r')

spreads2 = [sqrt(2/(fitresult0.s+1))*(180/pi), sqrt(2/(fitresult2.s+1))*(180/pi)]

%% When are cosine/sine products small?

om0 = 2*pi*0.1; om1 = 1.001*om0; om2 = 2*om0;

t = 0:0.0001:1e3;

n01 = cos(om0*t).*cos(om1*t);
n02 = cos(om0*t).*cos(om2*t);

recurrence_time = 1./(om1-om0)