%DISENO DE UNA RAT PARA EL DESPLIEGUE DEL TREN DE ATERRIZAJE AERONAVE MILITAR
%TRABAJO OBLIGATORIO 03 - SISTEMAS DE MOTOR, GIA, ETSIAE
%A. Cervigni, M. Piedrabuena, J. Rodriguez
syms t

%% SECCION 0 - DATOS DEL PROBLEMA Y PARAMETROS DEL ESTUDIO
%requisitos de diseno (datos del enunciado)
h = 1;          %m, distancia del CDG a la ariculación del tren
L = 2.22444;    %m, separacion entre puntos de rotación de tren y actuador
L_c = 1.22444;  %m, longitud de la camisa del actuador o carrera
L_p = 0.01;     %m, espesor del pistón
V_0 = 400/3.6;  %m/s, velocidad de vuelo
t_f = 10;       %s,  tiempo de extension de los trenes
D_p = 0.08;     %m, diámetro del pistón
D_r = 0.05;     %m, diámetro del empujador
C_dnlg = 2;     %Coef de resistencia aerodin. tren delantero
C_dmlg = 2.5;   %Coef de resistencia aerodin. tren principal
pr = 101330;    %Pa, presión en el lado de retracción
rho = 1.22;     %kg/m^3, densidad aire ambiente
S_nlg = 0.3;    %m^2, superficie frontal del tren delantero
S_mlg = 0.45;   %m^2, supeficie  frontal del tren principal
M_nlg = 60;     %kg, masa del tren frontal
M_mlg = 155;    %kg, masa del tren principal
rend_m = 0.97;  %rendimiento mecanico de la RAT
NACA = "23015"; %perfil aerodinamico de las palas de la RAT (teoria BEM)

ViscKin = 1.45*10^(-5);  %viscosidad cinematica del aire a nivel del mar (acorde densidad)
csound = 340.3;          %m/s, velocidad del sonido a nivel del mar (acorde densidad)
M_trans = 0.75;          %limitacion a Mach transonico en la punta del rotor (teoria 1D)

%parametros para el Xfoil
Re = 3*10^(5);  %numero de Reynolds medio para regimen a estudio
alpha0 = -20;    %angulo de ataque minimo de estudio del perfil
alphaf = 20;     %angulo de ataque maximo de estudio del perfil
dalpha = 0.25;   %variacion de angulo para generacion de polar
ITER = 1000;    %numero de iteraciones maximas del xfoil
filename = 'airfoildata.txt';

%datos para eleccion mediante teoria BEM (seccion 3)
B_s3 = 2;        %numero de palas

r0_s3 = 0.06;    %m, radio interior
rf_s3 = 0.18;    %m, radio exterior (punta)

N_int_s3 = 99;   %numero de intervalos de division de la pala

w0_s3 = 240;     %rad/s, velocidad angular mas pequena del estudio
wf_s3 = 300;     %rad/s, velocidad angular mas grande del estudio
Nw_s3 = 61;      %numero de divisiones intermedias entre las dos anteriores

theta0_s3 = -5;  %grados, angulo de paso mas pequeno de la optimizacion
thetaf_s3 = 70;  %grados, angulo de paso mas grande de la optimizacion
dtheta_s3 = 1;   %grados, aumento en cada iteracion del angulo de paso

%datos para estudio influencia parametros mediante teoria BEM (seccion 4)
B_s4 = 2;        %numero de palas

r0_s4 = 0.06;    %m, radio interior
rf0_s4 = 0.16;   %m, radio exterior mas pequeno del estudio
rff_s4 = 0.22;   %m, radio exterior mas grande del estudio
Nr_s4 = 13;      %numero de divisiones intermedias entre las dos anteriores  

N_int_s4 = 99;   %numero de intervalos de division de la pala

w0_s4 = 240;     %rad/s, velocidad angular mas pequena del estudio
wf_s4 = 300;     %rad/s, velocidad angular mas grande del estudio
Nw_s4 = 61;      %numero de divisiones intermedias entre las dos anteriores

theta0_s4 = -5;  %grados, angulo de paso mas pequeno de la optimizacion
thetaf_s4 = 70;  %grados, angulo de paso mas grande de la optimizacion
dtheta_s4 = 1;   %grados, aumento en cada iteracion del angulo de paso

%datos derivados
S_e = (pi*D_p^2)/4;        %m^2, area del embolo
S_r = (pi*D_r^2)/4;        %m^2, area del extensor
w_ext = (90/t_f)*(pi/180); %rad/s, velocidad angular de extension (cte)

% carga la polar (necesario tener el xfoil.exe en la misma carpeta o usar archivo externo)
XfoilAirfoilAnalysis(NACA, Re, alphaf, alpha0, dalpha,  ITER, filename);
[alfa_almacen, cl_almacen, cd_almacen] = XfoilDataRead(filename);

%% SECCION 1 - CINEMATICA Y DINAMICA
%cinematica del problema
gamma = w_ext*t;    
x = sqrt(h^2 + L^2 - 2*h*L*cos(gamma)) - L_c;
beta = asin((h*sin(gamma))/(L_c + x));
delta = asin((L*sin(gamma))/(L_c + x));

%cargas sobre el actuador del tren principal
W_mlg  = M_mlg*9.81;                       %peso del sistema aplicado al centro de masas
St_mlg = S_mlg*sin(gamma);                 %area frontal varible del sistema
D_mlg  = (1/2)*rho*St_mlg*V_0^(02)*C_dmlg; %resistencia aerodinamica del sistema

Fp_mlg = (D_mlg*sin(gamma) - W_mlg*cos(gamma))/...
  (cos(beta)*sin(gamma) + sin(beta)*cos(gamma));   %fuerza ejercida por el actuador hidraulico

Fpt_mlg = Fp_mlg*cos(beta)*sin(gamma);  %fuerza ejercida por el actuador en la linea de actuacion
root_mlg = double(vpasolve(Fp_mlg==0,t,[0 10]));

%presion hidraulica requerida por el sistema principal
ph_mlg = (Fp_mlg/S_e) + ((S_e-S_r)/S_e)*pr;

%cargas sobre el actuador del tren delantero
W_nlg  = M_nlg*9.81;                       %peso del sistema aplicado al centro de masas
St_nlg = S_nlg*sin(gamma);                 %area frontal varible del sistema
D_nlg  = (1/2)*rho*St_nlg*V_0^(02)*C_dnlg; %resistencia aerodinamica del sistema

Fp_nlg = (D_nlg*sin(gamma) - W_nlg*cos(gamma))/...
  (cos(beta)*sin(gamma) + sin(beta)*cos(gamma));   %fuerza ejercida por el actuador hidraulico

Fpt_nlg = Fp_nlg*cos(beta)*sin(gamma);  %fuerza ejercida por el actuador en la linea de actuacion
root_nlg = double(vpasolve(Fp_nlg==0,t,[0 10]));

%presion hidraulica requerida por el sistema delantero
ph_nlg = (Fp_nlg/S_e) + ((S_e-S_r)/S_e)*pr;

%potencia necesaria a suministrar por la RAT para el despliegue (teniendo en cuenta rendimiento mecanico)
x_dot = diff(x,t);
POW_NLG = Fp_nlg*x_dot;
POW_MLG = Fp_mlg*x_dot;

POW_MAX = double(subs((POW_NLG + 2* POW_MLG)/rend_m,10));

%% SECCION 2 - DISENO RAT MEDIANTE TEORIA 1D (con efecto de rotacion)
lambda = (M_trans*csound)/V_0; %valor de lambda para el estudio

%expresion para la potencia maxima obtenible dado delta (Glauert, 1935)
if (lambda >= 0.0) && (lambda <= 2.5)
    epsilon = -0.014*lambda^4 + 0.1433*lambda^3 - 0.5605*lambda^2 +...
        1.0502*lambda + 0.084;
else
    epsilon = 0.0003*lambda^3 - 0.008*lambda^2 + 0.0725*lambda + 0.763;
end

Cp_1D = (16*epsilon)/27;
D_1D = sqrt((8*POW_MAX)/(rho*(V_0)^3*pi*Cp_1D));
rpm_1D = ((2*lambda*V_0)/D_1D)*(30/pi);

%% SECCION 3 - DISENO RAT MEDIANTE TEORIA BEM
% valores para la optimizacion
dw_s3 = (wf_s3-w0_s3)/(Nw_s3-1);
Ntheta_s3 = ((thetaf_s3-theta0_s3)/dtheta_s3)+1;

%iniciacion variables almacen
Power_s3     = zeros(1,Nw_s3);
Cpsave_s3    = zeros(1,Nw_s3);
thetasave_s3 = zeros(1,Nw_s3);
Clsave_s3    = zeros(1,Nw_s3);

for n = 1:Nw_s3
    %parametros de la pala
    w  = w0_s3+dw_s3*(n-1);
    dr = (rf_s3-r0_s3)/N_int_s3;
    
    for p = 1:Ntheta_s3
        
        theta = (theta0_s3+dtheta_s3*(p-1))*((2*pi)/360);
        
        %inicializacion
        r = zeros(1, N_int_s3+1);
        c = zeros(1, N_int_s3+1);
        Pt = zeros(1, N_int_s3+1);
        Pn = zeros(1, N_int_s3+1);
        Cl = zeros(1, N_int_s3+1);
        
        %metodo BEM a la pala
        for i = 1:N_int_s3+1
            
            r(i) = r0_s3 + dr*(i-1);
            c(i) = 0.1 - (0.075/rf_s3)*(r(i));
            
            a_ax  = 0.0;
            a_tan = 0.0;
            a_ax_prev  = 10.0;
            a_tan_prev = 10.0;
            
            while (abs(a_ax-a_ax_prev)>0.01 && abs(a_tan-a_tan_prev)>0.01)
                
                a_ax_prev  = a_ax;
                a_tan_prev = a_tan;
                
                Ut        = w*r(i)*(1+a_tan_prev);
                Un        = V_0*(1-a_ax_prev);
                Vrel_norm = sqrt(Un.^2+Ut.^2);
                
                phi = atan2(Un,Ut);
                
                F = fTipLoss(B_s3,r(i),rf_s3,phi,true);
                
                alpha=(phi-theta)*(180/pi);
                
                sigma=(c(i)*B_s3*1.)/(2*pi*r(i));
                lambda_r = w*r(i)/V_0        ;
                
                if (alpha > alphaf || alpha<alpha0)
                    Cl = 0;
                    Cd = 0;
                else
                    Cl = interp1(alfa_almacen, cl_almacen, alpha, 'spline');
                    Cd = interp1(alfa_almacen, cd_almacen, alpha, 'spline');
                end
                
                cn=(Cl.*cos(phi)+Cd.*sin(phi));
                ct=(Cl.*sin(phi)-Cd.*cos(phi));
                
                Ct=Vrel_norm^2/V_0^2*sigma *cn;
                Cq=Vrel_norm^2/V_0^2*sigma *ct;
                
                [a_ax,a_tan] = fInductionCoefficients(a_ax_prev,Ct,Cq,F,lambda_r);
                
            end
            
            Pt(i)  =  0.5*rho*Vrel_norm.^2*c(i).*ct;
            
        end
        Power_prov = B_s3*trapz(r,r.*Pt)*w;
        
        if Power_prov > Power_s3(1,n)
            Power_s3(1,n)     = Power_prov
            Cpsave_s3(1,n)    = Power_prov/(0.5*rho*V_0^3*pi*(rf_s3)^2)
            Clsave_s3(1,n)    = Cl
            thetasave_s3(1,n) = theta*(180/pi)
        end
    end
end

%% SECCION 4 - ESTUDIO INFLUENCIA PARAMETROS MEDIANTE TEORIA BEM
% valores para la optimizacion
drit_s4   = (rff_s4-rf0_s4)/(Nr_s4-1);
dw_s4     = (wf_s4-w0_s4)/(Nw_s4-1);
Ntheta_s4 = ((thetaf_s4-theta0_s4)/dtheta_s4)+1;

%iniciacion variables almacen
Power_s4     = zeros(Nr_s4,Nw_s4);
Cpsave_s4    = zeros(Nr_s4,Nw_s4);
thetasave_s4 = zeros(Nr_s4,Nw_s4);
Clsave_s4    = zeros(Nr_s4,Nw_s4);

for m = 1:Nr_s4
    for n = 1:Nw_s4
        %parametros de la pala
        w  = w0_s4+dw_s4*(n-1);
        rf_s4 = rf0_s4+drit_s4*(m-1);
        dr = (rf_s4-r0_s4)/N_int_s4;
        
        for p = 1:Ntheta_s4
            
            theta = (theta0_s4+dtheta_s4*(p-1))*((2*pi)/360);

            %inicializacion
            r = zeros(1, N_int_s4+1);
            c = zeros(1, N_int_s4+1);
            Pt = zeros(1, N_int_s4+1);
            Pn = zeros(1, N_int_s4+1);
            Cl = zeros(1, N_int_s4+1);

            %metodo BEM a la pala
            for i = 1:N_int_s4+1

                r(i) = r0_s4 + dr*(i-1);
                c(i) = 0.1 - (0.075/rf_s4)*(r(i));

                a_ax  = 0.0;
                a_tan = 0.0;
                a_ax_prev  = 10.0;
                a_tan_prev = 10.0;

                while (abs(a_ax-a_ax_prev)>0.01 && abs(a_tan-a_tan_prev)>0.01)

                    a_ax_prev  = a_ax;
                    a_tan_prev = a_tan;

                    Ut        = w*r(i)*(1+a_tan_prev);
                    Un        = V_0*(1-a_ax_prev);
                    Vrel_norm = sqrt(Un.^2+Ut.^2);

                    phi = atan2(Un,Ut);

                    F = fTipLoss(B_s4,r(i),rf_s4,phi,true);

                    alpha=(phi-theta)*(180/pi);

                    sigma=(c(i)*B_s4*1.)/(2*pi*r(i));
                    lambda_r = w*r(i)/V_0        ;

                    if (alpha > (alphaf-0.25) || alpha<(alpha0+0.25))
                        Cl = 0;
                        Cd = 0;
                    else
                        Cl = interp1(alfa_almacen, cl_almacen, alpha, 'spline');
                        Cd = interp1(alfa_almacen, cd_almacen, alpha, 'spline');
                    end

                    cn=(Cl.*cos(phi)+Cd.*sin(phi));
                    ct=(Cl.*sin(phi)-Cd.*cos(phi));

                    Ct=Vrel_norm^2/V_0^2*sigma *cn;
                    Cq=Vrel_norm^2/V_0^2*sigma *ct;

                    [a_ax,a_tan] = fInductionCoefficients(a_ax_prev,Ct,Cq,F,lambda_r);

                end

                Pt(i)  =  0.5*rho*Vrel_norm.^2*c(i).*ct;

            end
            Power_prov = B_s4*trapz(r,r.*Pt)*w;
           
            if Power_prov > Power_s4(m,n)
                Power_s4(m,n)     = Power_prov
                Cpsave_s4(m,n)    = Power_prov/(0.5*rho*V_0^3*pi*(rf_s3)^2)
                Clsave_s4(m,n)    = Cl
                thetasave_s4(m,n) = theta*(180/pi)
            end
  
        end 
    end
end

%% SECCION 5 - REPRESENTACION GRAFICA
figure(1)
fplot(Fpt_mlg,[root_mlg t_f])

figure(2)
fplot(ph_mlg,[root_mlg t_f])

figure(3)
fplot(Fpt_nlg,[root_nlg t_f])

figure(4)
fplot(ph_nlg,[root_nlg t_f])

figure(5)
if root_mlg > root_nlg
    fplot(0.0, [0 root_nlg])
    hold on
    fplot(POW_NLG/rend_m, [root_nlg root_mlg])
    hold on
    fplot((POW_NLG + 2* POW_MLG)/rend_m, [root_mlg t_f])
    hold off
    grid on
elseif root_nlg > root_mlg
    fplot(0.0, [0 root_mlg])
    hold on
    fplot(2*POW_MLG/rend_m, [root_mlg root_nlg])
    hold on
    fplot((POW_NLG + 2* POW_MLG)/rend_m, [root_nlg t_f])
    hold off
    grid on
    
else
    fplot(0.0, [0 root_mlg])
    hold on
    fplot((POW_NLG + 2* POW_MLG)/rend_m, [root_mlg t_f])
    hold off
    grid on
end

figure(6)
plot(linspace(w0_s3,wf_s3,Nw_s3), Power_s4(1,:), 'blue')
hold on
plot(linspace(w0_s3,wf_s3,Nw_s3), Power_s4(5,:), 'blue')
hold on
plot(linspace(w0_s3,wf_s3,Nw_s3), Power_s4(9,:), 'blue')
hold on
plot(linspace(w0_s3,wf_s3,Nw_s3), Power_s4(13,:), 'blue')
hold off
grid on

figure(7)
stem3(linspace(w0_s4,wf_s4,Nw_s4), linspace(rf0_s4,rff_s4,Nr_s4), Power_s4, 'filled', 'LineStyle',':', 'Color','blue')

%% SECCION 6 - FUNCIONES DEL PROGRAMA
function [F] = fTipLoss(nB,r,R,phi,bTipLoss,varargin)
% - Compute tip-loss factor
    F=1;
    if bTipLoss && sin(phi)>0.01
        F=2/pi*acos(exp(-nB/2*(R-r)/(r*sin(phi)))); 
    end
end

function [a,aprime] = fInductionCoefficients(a_last,Ct,Cq,F,...
    lambda_r,varargin)
% - Compute a, aprime and the local thrust coefficient Ct
% - Perform High-thrust correction (e.g. a-Ct relation)
% - Perform relaxation on axial induction (only if steady simulation)
% - Perform wake-rotation correction
    [a,Ct] = fCorrectionHighThrust(Ct,F,varargin); % a-Ct relation 
    a = 0.3*a + (1-0.3)*a_last                   ; % Relaxation 
    aprime = Cq / (4*F*(1-a)*lambda_r)           ; % tangential induction
end

function [a,Ct] = fCorrectionHighThrust(Ct,F,varargin)
% - Returns a and Ct applying the High-Thrust correction
    k = [0.00,0.251163,0.0544955,0.0892074]; 
    Ctb = Ct./F;
    a = k(4)*Ctb.^3+k(3)*Ctb.^2+k(2)*Ctb+k(1);
end


function [alfa_almacen, cl_almacen, cd_almacen] = XfoilDataRead(filename)
% - Lee alpha, Cl y Cd de una polar generada con el programa Xfoil

    % abre el archivo de datos
    fid = fopen(filename, 'r');
    
    % lee el header, no se guarda (datos de partida)
    P = textscan(fid,' Calculated polar for: %[^\n]','Delimiter',' ','MultipleDelimsAsOne',true,'HeaderLines',3);
    P = textscan(fid, '%*s%*s%f%*s%f%s%s%s%s%s%s', 1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines', 2, 'ReturnOnError', false);
    P = textscan(fid, '%*s%*s%f%*s%*s%f%*s%f%*s%*s%f', 1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines', 0, 'ReturnOnError', false);

    % lee la informacion del archivo
    P = textscan(fid, '%f%f%f%f%f%f%f%*s%*s%*s%*s', 'Delimiter',  ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' , 4, 'ReturnOnError', false);
    
    % almacena la informacion
    alfa_almacen  = P{1}(:,1);
    cl_almacen    = P{2}(:,1);
    cd_almacen    = P{3}(:,1);
    
    %cierra el archivo temporal
    fclose(fid);
end

function XfoilAirfoilAnalysis(NACA, Re, alfaf, alfa0, dalfa, ITER, filename)
% - Genera la polar del perfil NACA dado el numero de Reynolds, el intervalo de
% angulos de ataque a estudiar, el numero de iteraciones y el propio codigo
% del perfil

    % borra el archivo de datos de la anterior sesion
    if (exist(filename,'file'))
        delete(filename);
    end
    
    % borra el archivo de configuracion del xfoil de la anterior sesion
    if (exist("xfoil_input.txt",'file'))
        delete("xfoil_input.txt");
    end
    
    % configura el Xfoil
    fid = fopen('xfoil_input.txt','w');   
    
    fprintf(fid,'NACA %s\n', NACA);
    fprintf(fid,'PLOP\n');
    fprintf(fid,'G F\n');
    fprintf(fid,'\n');
    fprintf(fid,'MDES\n');
    fprintf(fid,'FILT\n');
    fprintf(fid,'EXEC\n');
    fprintf(fid,'\n');
    fprintf(fid,'PANE\n');
    fprintf(fid,'OPER\n');
    fprintf(fid,'ITER %g\n', ITER);
    fprintf(fid,'RE %g\n', Re);
    fprintf(fid,'VISC %g\n', Re);
    fprintf(fid,'PACC\n');
    fprintf(fid,'%s\n', filename);
    fprintf(fid,'\n');
    fprintf(fid,'ASEQ\n');
    fprintf(fid,'%g\n', alfa0);
    fprintf(fid,'%g\n', alfaf);
    fprintf(fid,'%g\n', dalfa);
    fprintf(fid,'PACC\n');
    fprintf(fid,'VISC\n');
    fprintf(fid,'\n');
    fprintf(fid,'QUIT\n');

    % cierra el archivo de configuracion
    fclose(fid);
    
    % corre el Xfoil con los datos de entrada
    cmd = 'xfoil.exe < xfoil_input.txt';
    [status,result] = system(cmd);
end