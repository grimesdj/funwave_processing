function [Cm,Sm] = generate_wavemaker_Cm_Sm(Freq,theta,Hmo_each,mfreq,mtheta,fmin,fmax,h_gen,delta,fm,theta_input,sigma_theta_input,gamma_spec,DY,Nglob,alpha_c,PERIODIC,SMALL);

%% Create randomized phases:
    for kf=1:mfreq
        phi1(kf,1) = rand(1)*2*pi;
    end
    
%% First, we need the D_gen function...
%        CALL WK_NEW_EQUAL_DFREQ_IRREGULAR_WAVE &
%            (Nfreq,Ntheta,delta_WK,DEP_WK,FreqPeak,FreqMax,FreqMin,GammaTMA,&
%            Hmo,ThetaPeak,sigma_theta,rlamda_ir,beta_gen_ir,D_gen_ir,Phase_ir,&
%            Width_WK,omgn_ir,Periodic,DY,Nglob,Freq,alpha_c)
    ZERO=0.0;
    grav= 9.81;
        alpha=-0.39;
        alpha1=alpha+1.0/3.0;
        for kf = 1:mfreq
            ap = Hmo_each(1,kf)/sqrt(2.0)/2.0;
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

        
        rlamda(kf,1)=wkn*sin(theta(kf));
        beta_gen(kf)=80.0/delta^2/wave_length^2;
        rl_gen=wkn*cos(theta(kf));
        rI=sqrt(pi/beta_gen(kf))*exp(-rl_gen^2/4.0/beta_gen(kf));
        D_gen(kf,1)=2.0*ap*cos(theta(kf))...
            *(omgn(kf)^2-alpha1*grav*wkn^4*h_gen^3)...
            /(omgn(kf)*wkn*rI*(1.0-alpha*(wkn*h_gen)^2));

        end
% $$$         % WAVEMAKER WIDTH: This is unnecessary for our purposes...
% $$$         omgn_tmp=2.0*pi*fm;
% $$$         tb=omgn_tmp*omgn_tmp*h_gen/grav;
% $$$         tc=1.0+tb*alpha;
% $$$         wkn=sqrt((tc-sqrt(tc*tc-4.0*alpha1*tb))/(2.0*alpha1))/h_gen;
% $$$         width=delta*wave_length/2.0;
        
%% THEN:
%        CALL CALCULATE_NEW_Cm_Sm(Mloc,Nloc,DX,DY,Xc_WK,Ibeg,Jbeg,Nfreq,&
%            Ntheta,D_gen_ir,Phase_ir,Width_WK,rlamda_ir,beta_gen_ir,Cm,Sm)

    temp_index = 1;
    for kf = 2:mfreq;
        if Freq(kf)~=Freq(kf-1)
            temp_index = temp_index + 1;
        end
    end

    disp(['Number of distinct freqs: ',num2str(temp_index)])
    Cm = zeros(Nglob,mfreq);
    Sm = zeros(Nglob,mfreq);
    N  = Nglob;
    jjsta=1;
    Jbeg=1;
        for J=1:Nglob
            kkk = 0;
            kf = 1;
            %
            Cm(J,kf)=Cm(J,kf)...
                +D_gen(kf,1)...
                *cos(rlamda(kf,1)*((J-Jbeg+(jjsta-1))*DY-ZERO)+phi1(kf,1));
            %
            Sm(J,kf)=Sm(J,kf)...
                +D_gen(kf,1)...
                *sin(rlamda(kf,1)*((J-Jbeg+(jjsta-1))*DY-ZERO)+phi1(kf,1));
            %
            for kf=2:mfreq
                if Freq(kf)==Freq(kf-1)
                    kkk = kkk+1;
                    kf_temp = kf - kkk;
                    %
                    Cm(J,kf_temp)=Cm(J,kf_temp)...
                        +D_gen(kf,1)...
                        *cos(rlamda(kf,1)*((J-Jbeg+(jjsta-1))*DY-ZERO)+phi1(kf,1));
                        %
                    Sm(J,kf_temp)=Sm(J,kf_temp)...
                        +D_gen(kf,1)*...
                        sin(rlamda(kf,1)*((J-Jbeg+(jjsta-1))*DY-ZERO)+phi1(kf,1));
                else
                    kkk = 0;
                    %
                    Cm(J,kf)=Cm(J,kf)...
                        +D_gen(kf,1)*cos(rlamda(kf,1)*((J-Jbeg+(jjsta-1))*DY-ZERO)+phi1(kf,1));
                        %
                    Sm(J,kf)=Sm(J,kf)...
                        +D_gen(kf,1)*...
                        sin(rlamda(kf,1)*((J-Jbeg +(jjsta - 1))*DY-ZERO)+phi1(kf,1));
                end
            end
            end

end

