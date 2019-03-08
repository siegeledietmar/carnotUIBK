function [teta_s,teta_sM_WS,Rv,Rg,Rf,Rsif,U] = teta_ground_v1(t_ini,A,P,psi_slab,l_slab,psi_g,htop,z,VC,n,qw,zw,fall,lambdafi,dfi,lambdawi,dwi,lambdafci,dfci,radiator,teta_e,teta_FH_m,sim)

% Dteta_i = 0; z91 !!!
plotten = 0;

t0 = 8760*3600;
t = 0:3600:t0;

if nargin == 0
        plotten = 1;
        t_ini = 0;

        A = 7*10*2;
        P = (7+10)*2;
        psi_slab = 0.0;
        l_slab = P;
        psi_g = 0.0;
        htop = .5;
        z = 0; % Tiefe Keller
        VC = 0; % Volumen Keller
        n = 0.2;
        fall = 'slabonground';  % 'slabonground';
        lambdafi = [.15,1.4,0.04,1.33,.037,2.0];
        dfi = [0.13,0.08,0.04,.15,.22,.5];
        lambdawi = [.6,.7,.04,.7];
        dwi = [.015,.21,.2,.003];
        lambdafci = [];
        dfci = [];
        climate0 = load(['.\CARNOT_Strahour_final.dat']);
        teta_e = climate0(:,8);

        radiator = 1;
        teta_FH_m = 35;

end

        lambda = 2;
        rho = 2000;
        c = 1000;
        a = lambda/(rho*c);
        qw = 0.05; % Grundwasser
        zw = 5;
        G = 0.02; % 0.01 to 0.03 K/m ~ 0.06 W/m
        teta_i = 20;
        teta_grenz = 15; % 12 | 15
        teta_im = 20; % Winter
        teta_imS = 25; % Sommer

        M = 1:12;
        MD = [31,28,31,30,31,30,31,31,30,31,30,31];
        MDsum = [0 MD(1) sum(MD(1:2)) sum(MD(1:3)) sum(MD(1:4)) sum(MD(1:5)) sum(MD(1:6)) sum(MD(1:7)) sum(MD(1:8)) sum(MD(1:9)) sum(MD(1:10)) sum(MD(1:11)) sum(MD(1:12))];
        MDsum = MDsum*24;

        HGT = 0;
        HT = 0;

        for ii = 1:length(teta_e)
           if teta_e(ii) < teta_grenz
               HGT = HGT + (teta_i - teta_e(ii));
               HT = HT + 1;
           end
        end
        HGT = HGT / 24;
        HT = round(HT / 24);
        nHT = HT/365*12; 

        jj = 1;
        for j = 1: 365
            teta_e_d(j) = mean(teta_e(jj:jj+23));
            jj = jj+24;
        end
        for j = 1: length(M)
            teta_e_M(j) = mean(teta_e((MDsum(j)+1):(MDsum(j+1)+1)));
        end

%         teta_i_M = teta_im.*(teta_e_M<teta_grenz)+teta_imS.*(teta_e_M>=teta_grenz);
%         teta_iW = teta_im*ones(12,1);
%         teta_iS = teta_imS*ones(12,1);

    if radiator
        teta_i_M = teta_im*(teta_e_M<15)+teta_imS*(teta_e_M>=15);
        teta_iW = teta_im*ones(1,12);
        teta_iS = teta_imS*ones(1,12);
    else
        teta_i_M = teta_FH_m*(teta_e_M<15)+teta_imS*(teta_e_M>=15);
        teta_iW = teta_FH_m*ones(1,12);
        teta_iS = teta_imS*ones(1,12);
    end

%         Dteta_i = (max(teta_i_M) - min(teta_i_M))/2;
        Dteta_i = 0;
        teta_em = mean(teta_e_M);
        teta_em = teta_em + 1; % solar absorption
        Dteta_e = (max(teta_e_M) - min(teta_e_M))/2;
        B = A/(0.5 * P);
        lpsi_slab = sum(psi_slab .* l_slab);
% disp(['A = ' num2str(A) ' P = ' num2str(P) ' VC = ' num2str(VC) ' z = ' num2str(z) ' htop = ' num2str(htop) ' (l*psi)_slab = ' num2str(lpsi_slab) ' (P*psi)_g = ' num2str(psi_g*P)])

        if radiator
            Rsif = 0.17;
        else
            Rsif = 0.0;
        end
        switch fall
            case {'slabonground','slabongroundwithedgeinsulation'}
                Rsiw = 0.13;
                Rsew = 0.04;
                Rsef = 0.00;
            case {'heatedbasement'}
                Rsiw = 0.13;
                Rsew = 0.00; 
                Rsef = 0.00;
            case {'suspendedfloor','unheatedbasement'}
                Rsiw = 0.13;
                Rsew = 0.04;
                Rsef = 0.17;       
        end

        
        switch fall
            case 'heatedbasement'
                if sim
                    Rwi = dwi./lambdawi;
                    Rw = sum(Rwi);
                    Uw_tot = (Rsiw+Rsew+Rw)^-1;
                    Rwi(end) = [];
                    dwi(end) = [];
                    Rw = sum(Rwi);
                    Uw = (Rsiw+Rsew+Rw)^-1;
                    w = sum(dwi); % Wall thickness
                else
                    Rwi = dwi./lambdawi;
                    Rw = sum(Rwi);
                    Uw_tot = (Rsiw+Rsew+Rw)^-1;
                    Rw = sum(Rwi);
                    Uw = (Rsiw+Rsew+Rw)^-1;
                    w = sum(dwi); % Wall thickness  
                end
            otherwise
                Rwi = dwi./lambdawi;
                Rw = sum(Rwi);
                Uw = (Rsiw+Rsew+Rw)^-1;
                w = sum(dwi); % Wall thickness
        end
        
        if sim
            Rfi = dfi./lambdafi;
            Rf = sum(Rfi);
            Uf_tot = (Rsif + Rsef + Rf)^-1;
            Rfi(end) = [];
            dfi(end) = [];
            Rf = sum(Rfi);
            Uf = (Rsif + Rsef + Rf)^-1;
            f = sum(dfi); % Wall thickness
        else
            Rfi = dfi./lambdafi;
            Rf = sum(Rfi);
            Uf_tot = (Rsif + Rsef + Rf)^-1;
            Rf = sum(Rfi);
            Uf = (Rsif + Rsef + Rf)^-1;
            f = sum(dfi); % Wall thickness
        end

        Uff = Uf + lpsi_slab/A;
                
        Hreg = A*Uff;

        switch fall
        case {'unheatedbasement','heatedbasement'}
            Rfci = dfci./lambdafci;
            Rfc = sum(Rfci);
            Ufc = (Rsif + Rfc)^-1;
        otherwise
            Ufc = [];
        end

% disp(['U_w = ' num2str(Uw) ' W/(m^2 K) ' ' U_f = ' num2str(Uf) ' W/(m^2 K)' ' U_ff = ' num2str(Uff) ' W/(m^2 K)' ' Hreg = ' num2str(Hreg) ' W/K' ' U_fc = ' num2str(Ufc) ' W/(m^2 K)'])

        % dt = w + lambda*(Rsi + Rf + Rse);
        switch fall
        case ('unheatedbasement')
            dt = lambda/Ufc; % nach PHPP
        otherwise
            dt = lambda/Uff; % nach PHPP
            dw = [];
        end

        delta = (8760*3600/pi*a)^0.5;
        alpha = 1.5-12/(2*pi)*atan(dt/(dt+delta));
        beta = 1.5-0.42*log(delta/dt+1);

        Ubf = []; Ubw = []; Ux = []; Ug = [];
        % 'slabonground' | 'slabongroundwithedgeinsulation' | 'unheatedbasement' | 'suspendedfloor' | 'unheatedbasement'
        switch fall 
            case ('slabonground')
                insulation = 1;
                if insulation == 0
                    U = 2*lambda/(pi*B+dt)*ln(pi*B/dt+1);
                elseif insulation == 1
                    U = lambda/(0.457 * B + dt);
                elseif insulation == 2
                    R_g = 0.457 * B/lambda;
                    U = 1/(Rf+Rsif+Rsef+w/lambda+R_g);
                end
                Ubf = []; Ubw = []; 
                Hg = A*U + P*psi_g;
                Hpi = A*lambda/dt*(2/((1+delta/dt)^2+1))^0.5;
                Hpe = 0.37*P*lambda*log(delta/dt+1);
            case ('slabongroundwithedgeinsulation')
                horizontal = 1;
                insulation = 1;
                if insulation == 0
                    U = 2*lambda/(pi*B+dt)*ln(pi*B/dt+1);
                elseif insulation == 1
                    U = lambda/(0.457 * B + dt);
                elseif insulation == 2
                    R_g = 0.457 * B/lambda;
                    U = 1/(Rf+Rsif+Rsef+w/lambda+R_g);
                end
                Ubf = []; Ubw = [];
                D = 1;
                dn = 0.2;
                lambdan = 0.035;
                Rn = dn/lambdan;
                Rs = Rn - dn/lambda;
                ds = Rs *lambda;
                Hg = A*U + P*psi_g;
                Hpi = A*lambda/dt*(2/((1+delta/dt)^2+1))^0.5;
                if horizontal
                    psi_ge = -lambda/pi*((log(D/dt)+1)-log(D/(dt+ds)+1));
                    Hpe = 0.37 * P * lambda *((1-exp(-D/delta))*log(delta/(dt+ds)+1)+exp(-D/delta*log(delta/dt+1)));
                else
                    psi_ge = -lambda/pi*((log(2*D/dt)+1)-log(2*D/(dt+ds)+1));
                    Hpe = 0.37 * P * lambda *((1-exp(-2*D/delta))*log(delta/(dt+ds)+1)+exp(-2*D/delta*log(delta/dt+1)));
                end
            case ('suspendedfloor')
                h = 1; % Höhe der Bodenplatte über Boden bei Aufständerung
                dg = w+lambda*(Rsif+Rf+Rse);
                eps = 1; 
                v = 4; % Windgeschwindigkeit in m/s
                fw = 0.05;
                Ug = 2*lambda/(pi*B+dg)*log(pi*B/dg+1);
                Ux = 2*h*Uw/B+1450*eps*v*fw/B;
                U = (1/Uf + 1/(Ug+Ux))^-1;
                Ubf = []; Ubw = [];
                Hg = A*U + P*psi_g;
                Hpi = A*(1/Uf+1/(lambda/delta+Ux))^-1;
                Hpe = Uf * (0.37*P*lambda*log(delta/dg+1)+Ux*A)/(lambda/delta+Ux+Uf);
            case ('heatedbasement')
                if (dt+0.5*z)<B
%                     disp('(dt+0.5*z)<B');pause
                    Ubf = 2*lambda/(pi*B+dt+0.5*z)*log(pi*B/(dt+0.5*z)+1);
                else
%                     disp('(dt+0.5*z)>=B');pause
                    Ubf = lambda/(0.457*B+dt+0.5*z);
                end
                dw = lambda*(Rsiw+Rw+Rsew);
                Ubw = 2*lambda/(pi*z)*(1+0.5*dt/(dt+z))*log((z/dw+1));
                Uhf = (A*Ubf + z*P*Ubw)/(A+z*P);
                Hg = A*Ubf + z*P*Ubw+P*psi_g;
                Hpi = A*lambda/dt*(2/((1+delta/dt)^2+1))^0.5+z*P*lambda/dw*((1+delta/dw)^2+1)^0.5;
                Hpe = 0.37*P*lambda*(exp(-z/delta)*log(delta/dt+1)+2*(1-exp(-z/delta))*log(delta/dw+1));
                U = Uhf;
            case('unheatedbasement')
                dt = lambda/Ufc; % nach PHPP
                if (dt+0.5*z)<B
                    Ubf = 2*lambda/(pi*B+dt+0.5*z)*log(pi*B/(dt+0.5*z)+1);
                else
                    Ubf = lambda/(0.457*B+dt+0.5*z);
                end
                dw = lambda*(Rsiw+Rw+Rsew);
                Ubw = 2*lambda/(pi*z)*(1+0.5*dt/(dt+z))*log((z/dw+1));
                U = (1/Uf+A/(A*Ubf+z*P*Ubw+htop*P*Uw+0.33*n*VC))^(-1);
                Hg = A*U + P*psi_g;
                Hpi = (1/(A*Uf)+1/((A+z*P)*lambda/delta+htop*P*Uw+0.33*n*VC))^-1;
                Hpe = A*Uf*(0.37*P*lambda*(2-exp(-z/delta))*log(delta/dt+1)+htop*P*Uw+0.33*n*VC)/((A+z*P)*lambda/delta+htop*P*Uw+0.33*n*VC+A*Uf);
        end
        Hpe = Hpe + psi_g * P;
%         disp(['dt = ' num2str(dt) ' m, dw = ' num2str(dw) ' m'])

        rho_w_s = 1000;
        c_w_s = 4179;
        lc = lambda/(rho_w_s*c_w_s*qw/24/3600);
        G_w = 1 + 0.583*(1-2/pi*atan(23*lc/B))*exp(-3.75*zw/B)/(dt/B+0.35);

        Hg = Hg * G_w;
% disp(['U = ' num2str(U) ' W/(m^2 K), Ubf = ' num2str(Ubf) ' W/(m^2 K), Ubw = ' num2str(Ubw) ' W/(m^2 K), Hg = ' num2str(Hg) ' W/K, Hpe = ' num2str(Hpe) ' W/K, Hpi = ' num2str(Hpi) ' W/K'])
        Qdot_stat = Hg * (teta_im-teta_em);
%  disp(['Qdot_stat = ' num2str(Qdot_stat) ' W, Hg = ' num2str(Hg) ' W/K, teta_im = ' num2str(teta_im) ' °C, teta_em = ' num2str(teta_em) ' °C'])
        tau_e = 1;
        teta_eM = teta_em - Dteta_e*cos(2*pi*(M-tau_e)/12);

        tau_i = 7;
        teta_iM = teta_im + Dteta_i*cos(2*pi*(M-tau_i)/12);

        gamma = 12/(nHT*pi)*sin(nHT*pi/12)*cos(beta*pi/6);

        QdotM = Hg * (teta_im-teta_em) - ...
                Hpi * Dteta_i*cos(2*pi*(M-tau_i-alpha)/12) + ...
                Hpe * Dteta_e*cos(2*pi*(M-tau_e-beta)/12);

        QdotMW = Hg * (teta_iW-teta_em) - ...
                 Hpi * Dteta_i*cos(2*pi*(M-tau_i-alpha)/12) + ...
                 Hpe * Dteta_e*cos(2*pi*(M-tau_e-beta)/12);

        QdotMS = Hg * (teta_iS-teta_em) - ...
                 Hpi * Dteta_i*cos(2*pi*(M-tau_i-alpha)/12) + ...
                 Hpe * Dteta_e*cos(2*pi*(M-tau_e-beta)/12);

        Qdotm = Hg * (teta_im-teta_em) - ...
                Hpi * (teta_im-teta_iM) + ...
                Hpe * (teta_em-teta_eM);

         Qdot = Hg * (teta_im-teta_em)-gamma*Hpi*Dteta_i+gamma*Hpe*Dteta_e;
% disp('*')
% disp(['Qdot = ' num2str(Qdot) ' W, Hg = ' num2str(Hg) ' W/K, teta_im = ' num2str(teta_im) ' °C, teta_em = ' num2str(teta_em) ' °C,' ...
%       'gamma = ' num2str(gamma) ' W/K, Hpi = ' num2str(Hpi) ' W/K, Dteta_i = ' num2str(Dteta_i) ' K, Hpe = ' num2str(Hpe) ' W/K, Dteta_e = ' num2str(Dteta_e) ' K'])

         Qdot_harm = gamma*Hpe*Dteta_e;

         Q_HT = nHT * (Qdot_stat+Qdot_harm) * 8.76/12;

% disp(['Qdot_stat = ' num2str(Qdot_stat) ' W, Qdot_harm = ' num2str(Qdot_harm) ' W' ' Qdot = ' num2str(Qdot) ' W' ' Q_HT = ' num2str(Q_HT) ' kWh' ])

switch fall
    case ('heatedbasement')
%         UAeff = (Uff*A+psi_g*P+z*Uw*P);
%         UAeff = (Uf_tot*A+z*Uw_tot*P); % +psi_g*P
        if sim
            UAeff = Uf_tot*A+psi_g*P;
        else
            UAeff = Uf_tot*A+z*Uw_tot*P+psi_g*P;
        end
    otherwise
%         UAeff = (Uff*A+psi_g*P);
        UAeff = Uf_tot*A+psi_g*P;
end



         teta_sM = teta_im-(QdotM/UAeff);
         teta_sMW = teta_iW-(QdotMW/UAeff);
         teta_sMS = teta_iS-(QdotMS/UAeff);
disp(['UAeff = ' num2str(UAeff) ' W/K '])

teta_sM_WS = teta_sMW.*(teta_e_M<teta_grenz)+teta_sMS.*(teta_e_M>=teta_grenz);
Qdot_sM_WS = QdotMW.*(teta_e_M<teta_grenz)+QdotMS.*(teta_e_M>=teta_grenz);  
        teta_S = teta_em + G*z - Dteta_e * exp(-z*(pi/(a*t0))^0.5)*cos(z*(pi/(a*t0))^0.5-2*pi*t/t0);
        for j = 1: length(M)
            teta_S_M(j) = mean(teta_S((MDsum(j)+1):(MDsum(j+1)+1)));
        end


        disp(['A = ' num2str(A) ' P = ' num2str(P) ' VC = ' num2str(VC) ' z = ' num2str(z) ' htop = ' num2str(htop) ' (l*psi)_slab = ' num2str(lpsi_slab) ' (P*psi)_g = ' num2str(psi_g*P)])
        disp(['U_w = ' num2str(Uw) ' W/(m^2 K) ' ' U_f = ' num2str(Uf) ' W/(m^2 K)' ' U_ff = ' num2str(Uff) ' W/(m^2 K)' ' Hreg = ' num2str(Hreg) ' W/K' ' U_fc = ' num2str(Ufc) ' W/(m^2 K)'])
        disp(['dt = ' num2str(dt) ' m, delta = ' num2str(delta) ' m, alpha = ' num2str(alpha) ' m, beta = ' num2str(beta) ' M'])
        disp(['dt = ' num2str(dt) ' m, dw = ' num2str(dw) ' m'])
        disp(['U = ' num2str(U) ' W(m^2 K), Ubf = ' num2str(Ubf) ' W(m^2 K), Ubw = ' num2str(Ubw) ' W/(m^2 K), Hg = ' num2str(Hg) ' W/K, Hpe = ' num2str(Hpe) ' W/K, Hpi = ' num2str(Hpi) ' W/K'])
        disp(['Qdot_stat = ' num2str(Qdot_stat) ' W, Qdot_harm = ' num2str(Qdot_harm) ' W' ' Qdot = ' num2str(Qdot) ' W' ' Q_HT = ' num2str(Q_HT) ' kWh' ])


        a0 = [mean(teta_sM_WS),(max(teta_sM_WS)-min(teta_sM_WS))/2,0];
        a2 = fminsearch(@(a) norm((a0(1)+a0(2)*sin(M*2*pi/(12)+a))-teta_sM)',0); a1 = [mean(teta_sM_WS),(max(teta_sM_WS)-min(teta_sM_WS))/2,a2];
        disp([num2str(a1(1),'%4.1f') ' + (' num2str(a1(2),'%4.1f') ')*sin(2*pi/(8760*3600)*t+(' num2str(a1(3),'%4.1f') '))']);

        teta_s.bias = a1(1);
        teta_s.amp = a1(2);
        teta_s.freq = 2*pi*1/(8760*3600);
        % teta_s_phase = a1(3);
        teta_s.phase = a1(3)-2*pi/8760*t_ini;

        Rg = 0.5/lambda;
        Rv = max([0.01,1/U-Rsif-Rf-Rg]); disp(num2str(Rv))
%         Rv = max([0.01,1/Uf-Rsif-Rf-Rg]);
%         Rv = max([0.01,1/Uf_tot-Rsif-Rf-Rg]);
%         Rv = 1/Uf-Rsif-Rf-Rg;
%           Rv = 1/Uf-Rsif-Rf-Rg;
          Rv = 0.001;
          
        disp(['Uf = ' num2str(Uf) ' W/(m^2 K)' ' Uff = ' num2str(Uff) ' W/(m^2 K)'  ' Uf_tot = ' num2str(Uf_tot) ' W/(m^2 K)' ' U = ' num2str(U) ' W/(m^2 K)']) 
        
if plotten
figure(1);clf;subplot(2,1,1); hold on;
    plot([1:12],teta_sM_WS,'d-k')
    plot(t/3600/(8760/12),a1(1)+a1(2)*sin(2*pi/(8760*3600)*t+a1(3)),'k-');
    plot([1:12],teta_eM,'^--b')
    plot([1:12],teta_iM,'s-r')
    set(gca,'Xlim',[1,12])
    legend('soil','fit','e','i')
    grid on
    ylabel('\vartheta / [°C]','FontWeight','bold')
    subplot(2,1,2); hold on;
    plot([1:12],QdotM,'d-k')
    set(gca,'Xlim',[1,12])
    grid on
    ylabel('Q / [W]','FontWeight','bold')
end
