mfreq=1025;
mtheta=25;
fmin= 0.04;
fmax= 0.25;
h_gen=9;
delta=2;
fm = 0.1;
theta_input=0;
sigma_theta_input=10;
gamma_spec=3;
DY = 1;
Nglob=2999;
alpha_c=0;
PERIODIC=1;
SMALL=eps('single');


df = (fmax-fmin)/(mfreq-1.0);
    Ef = 0;
    for kf=1:mfreq
        Freq(kf) = fmin +(kf-1)*df;
    end

    % % --- directional wave
    sigma_theta=sigma_theta_input*pi/180.0;
    N_spec=20.0/sigma_theta;

    if mtheta==1 % 1D case
        theta(1) = theta_input*pi/180.0;
        AG(1) = 1.0;
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
            theta(kf) = (-1).^real(kf)*(-pi*1.0/3.0 + ...
                2.0/3.0*pi*(floor(real(ktheta_temp)/2.0 - ...
                0.5))/(real(mtheta)-1.0));   %Grimes narrowed +/- bounds

            theta(kf) = theta(kf) + theta_input*pi/180.0;   %new method
            if theta(kf)>0.5*pi
                theta(kf) = 0.5*pi;
            end
            if theta(kf)<-0.5*pi
                theta(kf) = -0.5*pi;
            end
            AG(kf) = 1.0/( 2.0*pi );
            
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