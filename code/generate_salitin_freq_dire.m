mfreq=1550;
mtheta=31;
fmin= 0.04;
fmax= 0.25;
h_gen=9;
delta=2;
fm = 0.1;
theta_input=0;
sigma_theta_input=02;
gamma_spec=3;
DY = 1;
Nglob=2999;
alpha_c=0;
PERIODIC=1;
SMALL=eps('single');

THETA_MAX_DEG = min(10*sigma_theta_input,60);
THETA_MAX = THETA_MAX_DEG*pi/180;

df = (fmax-fmin)/(mfreq-1.0);
    Ef = 0;
    for kf=1:mfreq
        Freq(kf) = fmin +(kf-1)*df;
    end

    % % --- directional wave
    sigma_theta=sigma_theta_input*pi/180.0;
    N_spec=20.0/sigma_theta;

    if mtheta==1 % 1D case
        theta(1:mfreq) = theta_input*pi/180.0;
        AG(1:mfreq) = 1.0;
    else
        [~,displace_theta] = min(abs(Freq-fm));
        idx_theta = mod(displace_theta(1),mtheta);
        for kf=1:mfreq %new method
            ktheta_temp = mod(kf-idx_theta,mtheta);
            if ktheta_temp<=0
                ktheta_temp = ktheta_temp + mtheta;
            end
            theta(kf) = (-1).^real(kf)*(-pi*1.0/2.0 + ...
                2.0/2.0*pi*(floor(real(ktheta_temp)/2.0 - ...
                0.5))/(real(mtheta)-1.0));   %new method
            theta(kf) = (-1).^real(kf)*(-THETA_MAX + ...
                2*THETA_MAX*(floor(real(ktheta_temp)/2.0 - ...
                0.5))/(real(mtheta)-1.0));   %Grimes narrowed +/- bounds

            theta(kf) = theta(kf) + theta_input*pi/180.0;   %new method
            if theta(kf)>THETA_MAX
                theta(kf) = THETA_MAX;
            end
            if theta(kf)<-THETA_MAX
                theta(kf) = -THETA_MAX;
            end
            AG(kf) = 1.0/( 2.0*pi );
% $$$             if abs(theta(kf))<pi/180 & abs(kf-displace_theta)<2
% $$$                 return
% $$$             end
            for k_n=1:N_spec
                AG(kf) = AG(kf)+ ...
                    (1.0/pi)*exp(-0.5*(real(k_n)*sigma_theta).^2) ...
  	                *cos(k_n*(theta(kf)-theta_input*pi/180.0));
            end
        end
        AG(:) = abs(AG(:));
    end

    theta0 = theta;

    % next correct the directions for periodic domain...
    alpha=-0.39;
    alpha1=alpha+1.0/3.0;
    grav=9.81;
    for kf = 1:mfreq
        omgn(kf)=2.0*pi*Freq(kf);
        tb=omgn(kf)*omgn(kf)*h_gen/grav;
        tc=1.0+tb*alpha;
        wkn=sqrt((tc-sqrt(tc*tc-4.0*alpha1*tb))/(2.0*alpha1))/h_gen;

        if wkn==0
            wkn=SMALL;
            C_phase=sqrt(grav*h_gen);
            wave_length=C_phase/fm;
        else
            C_phase=1.0/wkn*fm*2.0*pi;
            wave_length=C_phase/fm;
        end
        % for periodic boundary conditions
        theta_temp = theta(kf);
        if PERIODIC
            tmp1=wkn;
            if theta_temp>0
                tmp3=0;
                I=0;
                while (tmp3<theta_temp)
                    I=I+1;
                    tmp2=I*2.0*pi/DY/(Nglob-1.0);
                    if tmp2>=tmp1
                        theta_temp = theta_temp - 0.001;
                        if theta_temp<=0
                            theta_temp = 0.0;
                            break
                        end
                        break
                    else
                        % theta, based on rlamda=wkn*sin(theta)
                        tmp3=asin(tmp2/tmp1);
                    end
                end
                if tmp2<tmp1
                    tmp3=asin((I-1)*2.0*pi/DY/(Nglob-1.0)/tmp1);
                end
            elseif theta_temp<0
                tmp3=0;
                I=0;
                while (tmp3>theta_temp);
                    I=I+1;
                    tmp2=I*2.0*pi/DY/(Nglob-1.0);     % rlamda
                    if tmp2>=tmp1
                        theta_temp = theta_temp + 0.001;
                        if theta_temp>=0
                            theta_temp = 0.0;
                            break
                        end
                        break
                    else
                        % theta, based on rlamda=wkn*sin(theta)
                        tmp3=-asin(tmp2/tmp1);
                    end
                end
                if tmp2<tmp1
                    tmp3=-asin((I-1)*2.0*pi/DY/(Nglob-1.0)/tmp1);
                end
            elseif theta_temp==0.0
                tmp3 = theta_temp;
            end
        end
        theta(kf)=tmp3;
    end

    %% Wave heigh spectrum estimation:
    Hmo = 1.0;
    
        for kf=1:mfreq
        omiga_spec=2.0*pi*Freq(kf)*sqrt(h_gen/grav);
        phi=1.0-0.5*(2.0-omiga_spec).^2;
        if omiga_spec<=1.0, phi=0.5*omiga_spec.^2; end
        if omiga_spec>=2.0, phi=1.0; end
        sigma_spec=0.07;
        if Freq(kf)>fm, sigma_spec=0.09; end
        Etma(kf)=grav.^2*Freq(kf).^(-5)*(2.0*pi).^(-4)*phi...
            *exp(-5.0/4.0*(Freq(kf)/fm).^(-4))...
            *gamma_spec.^(exp(-(Freq(kf)/fm-1.0).^2/(2.0*sigma_spec.^2)));
        EnergyBin(kf) = Etma(kf)*df;
        Ef = Ef + EnergyBin(kf);
        end
        alpha_spec=Hmo.^2/16.0/Ef;
        correction_coeff = Ef/dot(AG,EnergyBin);
        for kf=1:mfreq
            AG(kf) = AG(kf) * correction_coeff;
            Hmo_each(1,kf)=4.0 *sqrt((alpha_spec*EnergyBin(kf)*AG(kf)));
        end


        figure, scatter(Freq,theta*180/pi,20,Hmo_each/2/sqrt(2),'filled')

% $$$         figure, scatter(Freq,theta*180/pi,20,log(Hmo_each),'filled')


        if mod(mfreq,mtheta)==0
            Hmo_array   = reshape(Hmo_each,mtheta,mfreq/mtheta);
            theta_array = reshape(theta0,mtheta,mfreq/mtheta);
            freq_array  = reshape(Freq,mtheta,mfreq/mtheta);
        end