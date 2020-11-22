%DISENO DE LOS CAMBIADORES DE CALOR DEL SISTEMA DE REFRIGERACION DE LOS ACEITES DE LA TURBINA LM2500 BASE DLE
%INTERCAMBIO DE CALOR ACEITE MINERAL - GAS NATURAL Y ACEITE SINTETICO - GAS NATURAL
%A. Cervigni, M. Piedrabuena, J. Rodriguez
syms T x

%% SECCION 0 - DATOS GLOBALES - REQUISITOS DE DISEÑO
R_u = 8.314; %J/molK

%datos de los componentes del gas natural
X_CH4 = 0.8425;    %fraccion molar metano
M_CH4 = 16;    %g/mol
Cp_CH4 = R_u*(4.568 + (-8.975*10^(-03))*T + (3.631*10^(-05))*T^2 + (-3.407*10^(-08))*T^3 + (1.091*10^(-11))*T^4);    %J/molK
k_CH4 = -0.00935+(1.4028*10^(-04))*T + (3.3180*10^(-08))*T^2;     %W/mK
nu_CH4 = (3.844+(4.0112*10^(-01))*T + (-1.4303*10^(-04))*T^2)*10^(-06);     %kg/ms

X_N2 = 0.0071;     %nitrogeno
M_N2 = 28;     %g/mol
Cp_N2 = R_u*(3.539 + (-0.261*10^(-03))*T + (0.007*10^(-05))*T^2 + (0.157*10^(-08))*T^3 + (-0.099*10^(-11))*T^4);     %J/molK
k_N2 = 0.00309 +(7.593*10^(-05))*T + (-1.1014*10^(-08))*T^2;      %W/mK
nu_N2 = (42.606+(4.7500*10^(-01))*T + (-9.8800*10^(-05))*T^2)*10^(-06);     %kg/ms

X_C2H6 = 0.1477;   %etano
M_C2H6 = 30;   %g/mol
Cp_C2H6 = R_u*(4.178 + (-4.427*10^(-03))*T + (5.660*10^(-05))*T^2 + (-6.651*10^(-08))*T^3 + (2.487*10^(-11))*T^4);   %J/molK
k_C2H6 = -0.01936 +(1.2547*10^(-04))*T + (3.8298*10^(-08))*T^2;   %W/mK
nu_C2H6 = (0.514+(3.3449*10^(-01))*T + (-7.1071*10^(-05))*T^2)*10^(-06);    %kg/ms

X_C3H8 = 0.0025;   %propano
M_C3H8 = 44;   %g/mol
Cp_C3H8 = R_u*(3.847 + (5.131^(-03))*T + (6.011*10^(-05))*T^2 + (-7.893*10^(-08))*T^3 + (3.079*10^(-11))*T^4);       %J/molK
k_C3H8 = -0.00869 +(6.6409*10^(-05))*T + (7.8760*10^(-08))*T^2;   %W/mK
nu_C3H8 = (-5.462+(3.2722*10^(-01))*T + (-1.0672*10^(-04))*T^2)*10^(-06);   %kg/ms

X_iC4H10 = 0.0001; %iso butano
M_iC4H10 = 58; %g/mol
Cp_iC4H10 = R_u*(3.351 + (17.883*10^(-03))*T + (5.477*10^(-05))*T^2 + (-8.099*10^(-08))*T^3 + (3.243*10^(-11))*T^4); %J/molK
k_iC4H10 = -0.00115 +(1.4943*10^(-05))*T + (1.4921*10^(-07))*T^2;  %W/mK
nu_iC4H10 = (-4.731+(2.9131*10^(-01))*T + (-8.0995*10^(-05))*T^2)*10^-(06); %kg/ms

X_nC4H10 = 0.0001; %normal butano
M_nC4H10 = 58; %g/mol
Cp_nC4H10 = R_u*(5.547 + (5.536*10^(-03))*T + (8.057*10^(-05))*T^2 + (-10.571*10^(-08))*T^3 + (4.134*10^(-11))*T^4); %J/molK
k_nC4H10 = -0.00182 +(1.9396*10^(-05))*T + (1.3818*10^(-07))*T^2;  %W/mK
nu_nC4H10 = (-4.946+(2.9001*10^(-01))*T + (-6.9665*10^(-05))*T^2)*10^(-06); %kg/ms

%datos del gas natural
vdot_g = 0.0612;                                           %m^3/s
pg_1 = 3400000;                                            %Pa
Tg_1 = 301.15;                                             %K
Tg_3_id = 353.15;                                          %K
M_g = X_CH4*M_CH4 + X_N2*M_N2 + X_C2H6*M_C2H6 + ...
    X_C3H8*M_C3H8 + X_iC4H10*M_iC4H10 + X_nC4H10*M_nC4H10; %g/mol
Cp_g = X_CH4*Cp_CH4 + X_N2*Cp_N2 + X_C2H6*Cp_C2H6 ...
    + X_C3H8*Cp_C3H8 + X_iC4H10*Cp_iC4H10 + X_nC4H10*Cp_nC4H10;
k_g = X_CH4*k_CH4 + X_N2*k_N2 + X_C2H6*k_C2H6 + ...
    X_C3H8*k_C3H8 + X_iC4H10*k_iC4H10 + X_nC4H10*k_nC4H10;
nu_g = X_CH4*nu_CH4 + X_N2*nu_N2 + X_C2H6*nu_C2H6 + ...
    X_C3H8*nu_C3H8 + X_iC4H10*nu_iC4H10 + X_nC4H10*nu_nC4H10;

% valores medios de las propiedades fluidas del gas natural
Cp_gm = double(int(Cp_g, T, Tg_1, Tg_3_id)/(Tg_3_id-Tg_1));   %J/molK
Cp_gm = (Cp_gm*1000)/M_g;                                     %J/KgK
k_gm = double(int(k_g, T, Tg_1, Tg_3_id)/(Tg_3_id-Tg_1));     %W/mK
nu_gm = double(int(nu_g, T, Tg_1, Tg_3_id)/(Tg_3_id-Tg_1));   %kg/ms
T_gm = 0.5*(Tg_1+Tg_3_id);                                    %K
rho_gm = (pg_1*M_g)/(R_u*1000*T_gm);                          %Kg/m^3
mdot_g = vdot_g*(pg_1/(((R_u*1000)/M_g)*T_gm));               %kg/s

%datos del aceite mineral (MLO)
Tmlo_1 = 340.15;                              %K
vdot_mlo = 0.00381884;                        %m^3/s
pmlo_1 = 5.4999474*10^5;                      %Pa
Cp_mlo = 3.6521739*T+809.3652174;             %J/(kg*K)
rho_mlo = -0.66666666666*T+1078.766666;       %kg/m^3
k_mlo = -7*10^(-05)*T + 0.1555;               %W/mK
v_mlo = 47313*10^(-06)*(T-273.15)^(-1.967);   %m^2/s
nu_mlo = rho_mlo*v_mlo;                       %Ns/m^2

%datos del aceite sintetico (TLO)
Tslo_1 = 366.48;                              %K
vdot_slo = 0.00381884;                        %m^3/s
pslo_1 = 5.4999474*10^5;                      %Pa
Cp_slo = 4.1821428*T+544;                     %J/(kg*K)
rho_slo = (-0.0007*T+1.077205)*1000;          %kg/m^3
k_slo = -7*10^(-05)*T + 0.134;                %W/mK
v_slo = 3086.8*10^(-06)*(T-273.15)^(-1.399);  %m^2/s
nu_slo = rho_slo*v_slo;                       %Ns/m^2

%% SECCION 1 - PREDIMENSIONADO
%proceso para el calculo de las temperaturas de los fluidos en las diferentes etapas en
%funcion de las razones de capacidad calorifica de cada fluido

%condiciones iniciales de la iteracion
Tmlo_2_id = Tmlo_1;
temp_mlo = Tmlo_2_id + 25;
Tslo_2_id = Tslo_1;
temp_slo = Tslo_2_id + 25;

%iteracion
while (abs(Tmlo_2_id-temp_mlo)>eps && abs(Tslo_2_id-temp_slo)>eps)

    temp_mlo = Tmlo_2_id; temp_slo = Tslo_2_id;

    rho_mlom = double( (subs(rho_mlo, T, Tmlo_2_id)+subs(rho_mlo, T, Tmlo_1)) / 2 ); %densidad lineal
    mdot_mlo = vdot_mlo*rho_mlom;

    rho_slom = double( (subs(rho_slo, T, Tslo_2_id)+subs(rho_slo, T, Tslo_1)) / 2 ); %densidad lineal
    mdot_slo = vdot_slo*rho_slom;

    Cp_mlom = double( (subs(Cp_mlo, T, Tmlo_2_id)+subs(Cp_mlo, T, Tmlo_1)) / 2 ); %Cp lineal
    Cp_slom = double( (subs(Cp_slo, T, Tslo_2_id)+subs(Cp_slo, T, Tslo_1)) / 2 ); %Cp lineal

    eqn = (Tg_3_id-x)/(x-Tg_1) - (((Tmlo_1-x)-((Tmlo_1 - (mdot_g*Cp_gm*(x-Tg_1))/(mdot_mlo*Cp_mlom)) - Tg_1))...
        / (log((Tmlo_1-x)/(Tmlo_1-(mdot_g*Cp_gm*(x-Tg_1))/(mdot_mlo*Cp_mlom)- Tg_1)))) /...
        (((Tslo_1-Tg_3_id)-((Tslo_1 - (mdot_g*Cp_gm*(Tg_3_id-x))/(mdot_slo*Cp_slom)) - x)) / ...
        (log((Tslo_1-Tg_3_id)/(Tslo_1-(mdot_g*Cp_gm*(Tg_3_id - x))/(mdot_slo*Cp_slom) - x))));

    Tg_2_id = double(vpasolve(eqn, x,327));

    Q1_id = mdot_g*Cp_gm*(Tg_2_id-Tg_1);
    Q2_id = mdot_g*Cp_gm*(Tg_3_id-Tg_2_id);
    Tmlo_2_id = Tmlo_1 - Q1_id/(mdot_mlo*Cp_mlom);
    Tslo_2_id = Tslo_1 - Q2_id/(mdot_slo*Cp_slom);
    varTm1_id = ((Tmlo_1-Tg_2_id)-(Tmlo_2_id-Tg_1))/log((Tmlo_1-Tg_2_id)/(Tmlo_2_id-Tg_1));
    varTm2_id = ((Tslo_1-Tg_3_id)-(Tslo_2_id-Tg_2_id))/log((Tslo_1-Tg_3_id)/(Tslo_2_id-Tg_2_id));

end

%calculo del area de transferencia de calor necesaria en cada cambiador en primera aproximacion
U0_1 = 200; %entre 200 y 400 para transf. gas a presion y liquido
U0_2 = 300;
U0_3 = 400;

F = 0.9;    %factor de correccion estimado para un shell and tube

A0_200 = Q1_id/(U0_1*F*varTm1_id);
A0_300 = Q1_id/(U0_2*F*varTm1_id);
A0_400 = Q1_id/(U0_3*F*varTm1_id);

%% SECCION 2 - ANALISIS DEL SISTEMA ELEGIDO
%debido a la carencia publica de datos especificos del cambiador optamos por estudiar el sistema a traves del metodo de Kern (simplificacion del
%Bell-Delaware) ya que precisa menos informacion geometrica para dar un resultado.

%calculo de las propiedades de los aceites en cada cambiador
k_mlom = double(int(k_mlo, T, Tmlo_1, Tmlo_2_id)/(Tmlo_2_id-Tmlo_1)); %W/mK
k_slom = double(int(k_slo, T, Tslo_1, Tslo_2_id)/(Tslo_2_id-Tslo_1)); %W/mK
nu_mlom = double(int(nu_mlo, T, Tmlo_1, Tmlo_2_id)/(Tmlo_2_id-Tmlo_1)); %kg/ms
nu_slom = double(int(nu_slo, T, Tslo_1, Tslo_2_id)/(Tslo_2_id-Tslo_1)); %kg/ms

%datos del cambiador de calor comunes a ambos a partir del catalogo
%configuracion "Cross Flow Baffle" con tubos disposicion triangular 30 grados
p_adm = 0.000145038*pg_1;                             %psi     %presion que va a soportar el interior del tubo
k_copper = 385;                                       %W/mK    %conductividad del material de la pared de los tubos
S_copper = 4830;                                      %psi     %allowable tensile stress
Y_copper = 0.4;                                       %        %wall thickness coefficient according ASME 31.3
E_copper = 1;                                         %        %quality factor for the piping according ASME 31.3
d_0 = 0.01905;                                        %m       %diametro externo de los tubos interiores
e_wall = 1.25*round((p_adm*d_0)/...
    (2*(S_copper*E_copper + p_adm*Y_copper)),4);      %m       %calculo del espesor de tubo segun norma ASME B31.3
d_i = d_0-2*e_wall;                                   %m       %diametro interno de los tubos interiores

%HX_data = (longitud, area intercambio calor)
HX_data = [0.6604, 0.648, 0.829, 45.3,  2; 0.6604, 0.648, 1.133, 55.7, 2; 0.6604, 0.673, 1.413, 66.1, 2;...
  0.6604, 0.705, 1.634, 76.5,  2; 0.6604, 0.705, 1.938, 86.9, 2; 0.6604, 0.724, 2.224, 97.3, 2;...
  0.6604, 0.724, 2.529, 107.7, 2;...
  0.6604, 0.578, 0.829, 44.2,  4; 0.6604, 0.578, 1.133, 54.3, 4; 0.6604, 0.603, 1.413, 64.5, 4;...
  0.6604, 0.635, 1.634, 74.6,  4; 0.6604, 0.635, 1.938, 84.7, 4; 0.6604, 0.654, 2.224, 94.9, 4;...
  0.6604, 0.654, 2.529, 105,   4;...
  0.6604, 0.546, 0.829, 41.5,  6; 0.6604, 0.546, 1.133, 51.0, 6; 0.6604, 0.572, 1.413, 60.5, 6;...
  0.6604, 0.603, 1.634, 70.0,  6; 0.6604, 0.603, 1.938, 79.5, 4; 0.6604, 0.622, 2.224, 89.1, 6;...
  0.6604, 0.622, 2.529, 98.6,  6;...
  0.7112, 0.651, 0.524, 42.4,  2; 0.7112, 0.676, 0.803, 54.5, 2; 0.7112, 0.708, 1.049, 66.6, 2;...
  0.7112, 0.708, 1.354, 78.8,  2; 0.7112, 0.727, 1.632, 90.9, 2; 0.7112, 0.727, 1.945, 103.1, 2;...
  0.7112, 0.727, 1.945, 115.2, 2;...
  0.7112, 0.613, 0.524, 41.5,  4; 0.7112, 0.638, 0.803, 53.4, 4; 0.7112, 0.670, 1.049, 65.3, 4;...
  0.7112, 0.670, 1.354, 77.2,  4; 0.7112, 0.689, 1.632, 89.1, 4; 0.7112, 0.689, 1.945, 101.1, 4;...
  0.7112, 0.721, 1.945, 113, 4;...
  0.7112, 0.549, 0.524, 39.0,  4; 0.7112, 0.575, 0.803, 50.2, 4; 0.7112, 0.606, 1.049, 61.4, 4;...
  0.7112, 0.606, 1.354, 72.6,  4; 0.7112, 0.625, 1.632, 83.7, 4; 0.7112, 0.625, 1.945, 94.9, 4;...
  0.7112, 0.657, 1.945, 106.1, 4;...
  0.7620, 0.683, 0.546, 50.1,  2; 0.7620, 0.740, 0.794, 64.1, 2; 0.762, 0.740, 1.065, 78,2;...
  0.7620, 0.759, 1.351, 92.1,  2; 0.7620, 0.759, 1.656, 105.8, 2; 0.762, 0.791, 1.929, 119.7, 2;...
  0.7620, 0.791, 2.234, 133.7, 2;...
  0.7620, 0.721, 0.546, 49.1,  4; 0.762, 0.676, 0.794, 62.8, 4; 0.7620, 0.676, 1.065, 76.4,  4;...
  0.7620, 0.695, 1.351, 90.1,  4; 0.762, 0.695, 1.656, 103.7, 4;0.7620, 0.727, 1.929, 117.4, 4;...
  0.7620, 0.727, 2.234, 131,   4;...
  0.7620, 0.587, 0.546, 46.4,  6; 0.7620, 0.638, 0.794, 59.4, 6; 0.7620, 0.638, 1.065, 72.3, 6;...
  0.7620, 0.664, 1.351, 85.2,  6; 0.7620, 0.664, 1.656, 98.1, 6; 0.7620, 0.689, 1.929, 111, 6;...
  0.7620, 0.689, 2.234, 123.8, 6];

Tg3_data = zeros(length(HX_data(:,1)), 1);
for i = 1:length(HX_data(:,1))
  for j = 1:length(HX_data(:,1))

    %datos del cambiador de calor primero a partir del catalogo
    D1 = HX_data(i,1);                  %m    %diametro interno de la carcasa
    L1 = HX_data(i,2)+ HX_data(i,3);    %m    %longitud del tubo
    E1 = HX_data(i,3);                  %m    %total baffle length
    A1 = HX_data(i,4);                  %m^2  %area de transferencia de calor
    Np1 = HX_data(i,5);                 %numero de tubos que pasan

    %datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
    %(datos estimados minimos para resolver por Kern)
    Pt1 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
    B1 = 0.4*D1;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
    C1 = Pt1-d_0;                     %espaciado libre entre tubos
    Nt1 = round(A1/(pi*d_0*L1), 0);   %numero de tubos interiores
    Nb1 = round((E1/B1)-1, 0);        %Numero de deflectores en el cambiador

    %datos del cambiador de calor primero a partir del catalogo
    D2 = HX_data(j,1);                  %m    %diametro interno de la carcasa
    L2 = HX_data(j,2)+ HX_data(j,3);    %m    %longitud del tubo
    E2 = HX_data(j,3);                  %m    %total baffle length
    A2 = HX_data(j,4);                  %m^2  %area de transferencia de calor
    Np2 = HX_data(j,5);                 %numero de tubos que pasan

    %datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
    %(datos estimados minimos para resolver por Kern)
    Pt2 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
    B2 = 0.4*D2;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
    C2 = Pt2-d_0;                     %espaciado libre entre tubos
    Nt2 = round(A2/(pi*d_0*L2), 0);   %numero de tubos interiores
    Nb2 = round((E2/B2)-1, 0);        %Numero de deflectores en el cambiador

    %desarrollo del metodo de KERN para el primer cambiador de calor
    D_e1  = 4*(((Pt1)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
        4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
    A_s1  = (D1*C1*B1)/(Pt1);                     %bundle crossflow area
    G_s1  = mdot_mlo/A_s1;                        %shell side mass velocity

    % shell-side
    T_w1=1/2*((Tg_1+Tg_2_id)/2+(Tmlo_1+Tmlo_2_id)/2);    %temperatura estimada en la pared de tubo
    nu_mlow = double(subs(nu_mlo,T,T_w1));               %viscosidad dinamica a la temperatura anterior
    Re_s1 = (G_s1*D_e1)/nu_mlom;                         %numero de reynolds del aceite en el shell
    phi_s1 = (nu_mlom/nu_mlow)^(0.14);
    h_mlo=(k_mlom/D_e1)*0.36*(D_e1*G_s1/nu_mlom)^(0.55)*...
        (Cp_mlom*nu_mlom/k_mlom)^(1/3)*phi_s1;           %coeficiente de transferencia de calor por convección del aceite

    % tube-side
    A_tp1=(pi*d_i^(02)/4)*(Nt1/Np1);                       %area total de entrada del gas a los tubos
    u_g1=(mdot_g)/(rho_gm*A_tp1);                          %velocidad del gas en los tubos
    Re_g1=(rho_gm*u_g1*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
    Pr_g1=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
    Nu_g1=(((1.58*log(Re_g1)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
        (Re_g1-1000)*Pr_g1)/(1+12.7*((1.58*...             %para flujo turbulento
        log(Re_g1)-3.28)^(-02)/2)^(0.5)*(Pr_g1^(2/3)-1));
    h_g1=Nu_g1*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas

    % resolucion
    U_1 = 1/((d_0/(d_i*h_g1))+(d_0*log(d_0/d_i)/(2*k_copper))+(1/h_mlo)); %coeficiente total de transferencia de calor
    eqn2 = (A1*U_1*(((Tmlo_1-Tg_1)-(Tmlo_2_id-x))/(log((Tmlo_1-Tg_1)/...
        (Tmlo_2_id-x))))/(mdot_g*Cp_gm)+Tg_1-x);
    Tg_2 = double(vpasolve(eqn2,x,353));                                  %temperatura del gas tras el primer cambiador

    %desarrollo del metodo de KERN para el segundo cambiador de calor
    D_e2  = 4*(((Pt2)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
        4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
    A_s2  = (D2*C2*B2)/(Pt2);                     %bundle crossflow area
    G_s2  = mdot_slo/A_s2;                        %shell side mass velocity

    % shell-side
    T_w2 = 1/2*((Tg_2_id+Tg_3_id)/2+(Tslo_1+Tslo_2_id)/2);  %temperatura estimada en la pared de tubo
    nu_slow = double(subs(nu_slo,T,T_w2));                  %viscosidad dinamica a la temperatura anterior
    Re_s2 = (G_s2*D_e2)/nu_slom;                            %numero de reynolds del aceite en el shell
    phi_s2 = (nu_slom/nu_slow)^(0.14);
    h_slo=(k_slom/D_e2)*0.36*(D_e2*G_s2/nu_slom)^(0.55)*...
        (Cp_slom*nu_slom/k_slom)^(1/3)*phi_s2;           %coeficiente de transferencia de calor por convección del aceite

    % tube-side
    A_tp2=(pi*d_i^(02)/4)*(Nt2/Np2);                       %area total de entrada del gas a los tubos
    u_g2=(mdot_g)/(rho_gm*A_tp2);                          %velocidad del gas en los tubos
    Re_g2=(rho_gm*u_g2*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
    Pr_g2=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
    Nu_g2=(((1.58*log(Re_g2)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
        (Re_g2-1000)*Pr_g2)/(1+12.7*((1.58*...             %para flujo turbulento
        log(Re_g2)-3.28)^(-02)/2)^(0.5)*(Pr_g2^(2/3)-1));
    h_g2=Nu_g2*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas                                       
    
    % resolucion
    U_2=1/((d_0/(d_i*h_g2))+(d_0*log(d_0/d_i)/(2*k_copper))+(1/h_slo));  %coeficiente total de transferencia de calor                                          
    eqn3 = (A2*U_2*(((Tslo_1-Tg_2)-(Tslo_2_id-x))/(log((Tslo_1-Tg_2)/...
        (Tslo_2_id-x))))/(mdot_g*Cp_gm)+Tg_2-x);
    Tg3_data(i,j) = double(vpasolve(eqn3,x,353))                        %temperatura del gas tras el segundo cambiador

  end
end

%% SECCION 3 - RESULTADO FINAL DEL CAMBIADOR LIMPIO
%datos del cambiador de calor primero a partir del catalogo (coger fila y
%columna de Tg3_data que de la temperatura final deseada, fila(i) para
%cambiador 1 y columna(j) para el cambiador 2)
D1 = HX_data(33,1);                  %m    %diametro interno de la carcasa
L1 = HX_data(33,2)+ HX_data(33,3);   %m    %longitud del tubo
E1 = HX_data(33,3);                  %m    %total baffle length
A1 = HX_data(33,4);                  %m^2  %area de transferencia de calor
Np1 = HX_data(33,5);                 %numero de tubos que pasan

%datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
%(datos estimados minimos para resolver por Kern)
Pt1 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
B1 = 0.4*D1;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
C1 = Pt1-d_0;                     %espaciado libre entre tubos
Nt1 = round(A1/(pi*d_0*L1), 0);   %numero de tubos interiores
Nb1 = round((E1/B1)-1, 0);        %Numero de deflectores en el cambiador

%datos del cambiador de calor primero a partir del catalogo
D2 = HX_data(56,1);                  %m    %diametro interno de la carcasa
L2 = HX_data(56,2)+ HX_data(56,3);   %m    %longitud del tubo
E2 = HX_data(56,3);                  %m    %total baffle length
A2 = HX_data(56,4);                  %m^2  %area de transferencia de calor
Np2 = HX_data(56,5);                 %numero de tubos que pasan

%datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
%(datos estimados minimos para resolver por Kern)
Pt2 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
B2 = 0.4*D2;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
C2 = Pt2-d_0;                     %espaciado libre entre tubos
Nt2 = round(A2/(pi*d_0*L2), 0);   %numero de tubos interiores
Nb2 = round((E2/B2)-1, 0);        %Numero de deflectores en el cambiador

%desarrollo del metodo de KERN para el primer cambiador de calor
D_e1  = 4*(((Pt1)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
    4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
A_s1  = (D1*C1*B1)/(Pt1);                     %bundle crossflow area
G_s1  = mdot_mlo/A_s1;                        %shell side mass velocity

% shell-side
T_w1=1/2*((Tg_1+Tg_2_id)/2+(Tmlo_1+Tmlo_2_id)/2);    %temperatura estimada en la pared de tubo
nu_mlow = double(subs(nu_mlo,T,T_w1));               %viscosidad dinamica a la temperatura anterior
Re_s1 = (G_s1*D_e1)/nu_mlom;                         %numero de reynolds del aceite en el shell
phi_s1 = (nu_mlom/nu_mlow)^(0.14);
h_mlo=(k_mlom/D_e1)*0.36*(D_e1*G_s1/nu_mlom)^(0.55)*...
    (Cp_mlom*nu_mlom/k_mlom)^(1/3)*phi_s1;           %coeficiente de transferencia de calor por convección del aceite

% tube-side
A_tp1=(pi*d_i^(02)/4)*(Nt1/Np1);                       %area total de entrada del gas a los tubos
u_g1=(mdot_g)/(rho_gm*A_tp1);                          %velocidad del gas en los tubos
Re_g1=(rho_gm*u_g1*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
Pr_g1=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
Nu_g1=(((1.58*log(Re_g1)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
    (Re_g1-1000)*Pr_g1)/(1+12.7*((1.58*...             %para flujo turbulento
    log(Re_g1)-3.28)^(-02)/2)^(0.5)*(Pr_g1^(2/3)-1));
h_g1=Nu_g1*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas

% resolucion
U_1 = 1/((d_0/(d_i*h_g1))+(d_0*log(d_0/d_i)/(2*k_copper))+(1/h_mlo)); %coeficiente total de transferencia de calor
eqn2 = (A1*U_1*(((Tmlo_1-Tg_1)-(Tmlo_2_id-x))/(log((Tmlo_1-Tg_1)/...
    (Tmlo_2_id-x))))/(mdot_g*Cp_gm)+Tg_1-x);
Tg_2 = double(vpasolve(eqn2,x,353));                                  %temperatura del gas tras el primer cambiador

%desarrollo del metodo de KERN para el segundo cambiador de calor
D_e2  = 4*(((Pt2)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
    4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
A_s2  = (D2*C2*B2)/(Pt2);                     %bundle crossflow area
G_s2  = mdot_slo/A_s2;                        %shell side mass velocity

% shell-side
T_w2 = 1/2*((Tg_2_id+Tg_3_id)/2+(Tslo_1+Tslo_2_id)/2);  %temperatura estimada en la pared de tubo
nu_slow = double(subs(nu_slo,T,T_w2));                  %viscosidad dinamica a la temperatura anterior
Re_s2 = (G_s2*D_e2)/nu_slom;                            %numero de reynolds del aceite en el shell
phi_s2 = (nu_slom/nu_slow)^(0.14);
h_slo=(k_slom/D_e2)*0.36*(D_e2*G_s2/nu_slom)^(0.55)*...
    (Cp_slom*nu_slom/k_slom)^(1/3)*phi_s2;           %coeficiente de transferencia de calor por convección del aceite

% tube-side
A_tp2=(pi*d_i^(02)/4)*(Nt2/Np2);                       %area total de entrada del gas a los tubos
u_g2=(mdot_g)/(rho_gm*A_tp2);                          %velocidad del gas en los tubos
Re_g2=(rho_gm*u_g2*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
Pr_g2=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
Nu_g2=(((1.58*log(Re_g2)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
    (Re_g2-1000)*Pr_g2)/(1+12.7*((1.58*...             %para flujo turbulento
    log(Re_g2)-3.28)^(-02)/2)^(0.5)*(Pr_g2^(2/3)-1));
h_g2=Nu_g2*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas                                       

% resolucion
U_2=1/((d_0/(d_i*h_g2))+(d_0*log(d_0/d_i)/(2*k_copper))+(1/h_slo));  %coeficiente total de transferencia de calor                                          
eqn3 = (A2*U_2*(((Tslo_1-Tg_2)-(Tslo_2_id-x))/(log((Tslo_1-Tg_2)/...
    (Tslo_2_id-x))))/(mdot_g*Cp_gm)+Tg_2-x);
Tg_3 = double(vpasolve(eqn3,x,353));                                 %temperatura del gas tras el segundo cambiador

% calculo del resto de propiedades de los fluidos
Q1 = mdot_g*Cp_gm*(Tg_2-Tg_1);
Q2 = mdot_g*Cp_gm*(Tg_3-Tg_2);
Tmlo_2 = Tmlo_1 - Q1/(mdot_mlo*Cp_mlom);
Tslo_2 = Tslo_1 - Q2/(mdot_slo*Cp_slom);
varTm1 = ((Tmlo_1-Tg_2)-(Tmlo_2-Tg_1))/log((Tmlo_1-Tg_2)/(Tmlo_2-Tg_1));
varTm2 = ((Tslo_1-Tg_3)-(Tslo_2-Tg_2))/log((Tslo_1-Tg_3)/(Tslo_2-Tg_2));
d_Pt1=(4*(1.58*log(Re_g1)-3.28)^(-02)*(L1*Np1)/(d_i) + 4*Np1)*(rho_gm*u_g1^(02))/2;
d_Pt2=(4*(1.58*log(Re_g2)-3.28)^(-02)*(L2*Np2)/(d_i) + 4*Np2)*(rho_gm*u_g2^(02))/2;
d_Ps1=((exp(0.576-0.19*log(Re_s1)))*G_s1^(02)*(Nb1+1)*D1)/(2*rho_mlom*D_e1*phi_s1);
d_Ps2=((exp(0.576-0.19*log(Re_s2)))*G_s2^(02)*(Nb2+1)*D2)/(2*rho_slom*D_e2*phi_s2);

%% SECCION 4 - RESULTADO FINAL DEL CAMBIADOR CON INCRUSTACION
%factor de incrustacion para el analisis de la suciedad en el cambiador
Rf_i = 0.000352;       %m^2*K/W   %factor de incrustacion en la pared interior del tubo
Rf_0 = 0.000176;       %m^2*K/W   %factor de incrustacion en la carcasa

%datos del cambiador de calor primero a partir del catalogo (coger fila y
%columna de Tg3_data que de la temperatura final deseada, fila(i) para
%cambiador 1 y columna(j) para el cambiador 2)
D1 = HX_data(33,1);                  %m    %diametro interno de la carcasa
L1 = HX_data(33,2)+ HX_data(33,3);   %m    %longitud del tubo
E1 = HX_data(33,3);                  %m    %total baffle length
A1 = HX_data(33,4);                  %m^2  %area de transferencia de calor
Np1 = HX_data(33,5);                 %numero de tubos que pasan

%datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
%(datos estimados minimos para resolver por Kern)
Pt1 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
B1 = 0.4*D1;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
C1 = Pt1-d_0;                     %espaciado libre entre tubos
Nt1 = round(A1/(pi*d_0*L1), 0);   %numero de tubos interiores
Nb1 = round((E1/B1)-1, 0);        %Numero de deflectores en el cambiador

%datos del cambiador de calor primero a partir del catalogo
D2 = HX_data(56,1);                  %m    %diametro interno de la carcasa
L2 = HX_data(56,2)+ HX_data(56,3);   %m    %longitud del tubo
E2 = HX_data(56,3);                  %m    %total baffle length
A2 = HX_data(56,4);                  %m^2  %area de transferencia de calor
Np2 = HX_data(56,5);                 %numero de tubos que pasan

%datos del cambiador de calor derivados de los anteriores a traves de la norma ASME
%(datos estimados minimos para resolver por Kern)
Pt2 = 1.25*d_0;                   %valor tipico, habitual entre 1.25 y 1.5
B2 = 0.4*D2;                      %espaciado entre deflectores - valor tipico, habitual entre 0.4 y 0.6
C2 = Pt2-d_0;                     %espaciado libre entre tubos
Nt2 = round(A2/(pi*d_0*L2), 0);   %numero de tubos interiores
Nb2 = round((E2/B2)-1, 0);        %Numero de deflectores en el cambiador

%desarrollo del metodo de KERN para el primer cambiador de calor
D_e1  = 4*(((Pt1)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
    4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
A_s1  = (D1*C1*B1)/(Pt1);                     %bundle crossflow area
G_s1  = mdot_mlo/A_s1;                        %shell side mass velocity

% shell-side
T_w1=1/2*((Tg_1+Tg_2_id)/2+(Tmlo_1+Tmlo_2_id)/2);    %temperatura estimada en la pared de tubo
nu_mlow = double(subs(nu_mlo,T,T_w1));               %viscosidad dinamica a la temperatura anterior
Re_s1 = (G_s1*D_e1)/nu_mlom;                         %numero de reynolds del aceite en el shell
phi_s1 = (nu_mlom/nu_mlow)^(0.14);
h_mlo=(k_mlom/D_e1)*0.36*(D_e1*G_s1/nu_mlom)^(0.55)*...
    (Cp_mlom*nu_mlom/k_mlom)^(1/3)*phi_s1;           %coeficiente de transferencia de calor por convección del aceite

% tube-side
A_tp1=(pi*d_i^(02)/4)*(Nt1/Np1);                       %area total de entrada del gas a los tubos
u_g1=(mdot_g)/(rho_gm*A_tp1);                          %velocidad del gas en los tubos
Re_g1=(rho_gm*u_g1*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
Pr_g1=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
Nu_g1=(((1.58*log(Re_g1)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
    (Re_g1-1000)*Pr_g1)/(1+12.7*((1.58*...             %para flujo turbulento
    log(Re_g1)-3.28)^(-02)/2)^(0.5)*(Pr_g1^(2/3)-1));
h_g1=Nu_g1*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas

% resolucion
U_1s = 1/((d_0/(d_i*h_g1))+((d_0*Rf_i)/d_i)+...
    (d_0*log(d_0/d_i)/(2*k_copper))+Rf_0+(1/h_mlo));                   %coeficiente total de transferencia de calor
eqn4 = (A1*U_1s*(((Tmlo_1-Tg_1)-(Tmlo_2_id-x))/(log((Tmlo_1-Tg_1)/...
    (Tmlo_2_id-x))))/(mdot_g*Cp_gm)+Tg_1-x);
Tg_2s = double(vpasolve(eqn4,x,353));                                  %temperatura del gas tras el primer cambiador

%desarrollo del metodo de KERN para el segundo cambiador de calor
D_e2  = 4*(((Pt2)^(02)*sqrt(3)/...            %Diametro equivalente para una distribución triangular de
    4 - pi*(d_0)^(02)/8)/(pi*(d_0)/2));       %tubos (4*free-flow area)/wetted perimeter
A_s2  = (D2*C2*B2)/(Pt2);                     %bundle crossflow area
G_s2  = mdot_slo/A_s2;                        %shell side mass velocity

% shell-side
T_w2 = 1/2*((Tg_2_id+Tg_3_id)/2+(Tslo_1+Tslo_2_id)/2);  %temperatura estimada en la pared de tubo
nu_slow = double(subs(nu_slo,T,T_w2));                  %viscosidad dinamica a la temperatura anterior
Re_s2 = (G_s2*D_e2)/nu_slom;                            %numero de reynolds del aceite en el shell
phi_s2 = (nu_slom/nu_slow)^(0.14);
h_slo=(k_slom/D_e2)*0.36*(D_e2*G_s2/nu_slom)^(0.55)*...
    (Cp_slom*nu_slom/k_slom)^(1/3)*phi_s2;           %coeficiente de transferencia de calor por convección del aceite

% tube-side
A_tp2=(pi*d_i^(02)/4)*(Nt2/Np2);                       %area total de entrada del gas a los tubos
u_g2=(mdot_g)/(rho_gm*A_tp2);                          %velocidad del gas en los tubos
Re_g2=(rho_gm*u_g2*d_i)/(nu_gm);                       %numero de Reynolds del gas en los tubos
Pr_g2=(Cp_gm*nu_gm)/(k_gm);                            %numero de Prandt del gas en los tubos
Nu_g2=(((1.58*log(Re_g2)-3.28)^(-02)/2)*...            %numero de Nusselt caldculado por la correlación de Gnielinski
    (Re_g2-1000)*Pr_g2)/(1+12.7*((1.58*...             %para flujo turbulento
    log(Re_g2)-3.28)^(-02)/2)^(0.5)*(Pr_g2^(2/3)-1));
h_g2=Nu_g2*k_gm/d_i;                                   %coeficiente de transferencia de calor por convección del gas                                      

% resolucion
U_2s=1/((d_0/(d_i*h_g2))+((d_0*Rf_i)/d_i)+...
    (d_0*log(d_0/d_i)/(2*k_copper))+Rf_0+(1/h_slo));                  %coeficiente total de transferencia de calor            
eqn5 = (A2*U_2s*(((Tslo_1-Tg_2s)-(Tslo_2_id-x))/(log((Tslo_1-Tg_2s)/...
    (Tslo_2_id-x))))/(mdot_g*Cp_gm)+Tg_2s-x);
Tg_3s = double(vpasolve(eqn5,x,353));                                 %temperatura del gas tras el segundo cambiador

% calculo del resto de propiedades de los fluidos
Q1s = mdot_g*Cp_gm*(Tg_2s-Tg_1);
Q2s = mdot_g*Cp_gm*(Tg_3s-Tg_2s);
Tmlo_2s = Tmlo_1 - Q1s/(mdot_mlo*Cp_mlom);
Tslo_2s = Tslo_1 - Q2s/(mdot_slo*Cp_slom);
varTm1s = ((Tmlo_1-Tg_2s)-(Tmlo_2s-Tg_1))/log((Tmlo_1-Tg_2s)/(Tmlo_2s-Tg_1));
varTm2s = ((Tslo_1-Tg_3s)-(Tslo_2s-Tg_2s))/log((Tslo_1-Tg_3s)/(Tslo_2s-Tg_2));
d_Pt1s=(4*(1.58*log(Re_g1)-3.28)^(-02)*(L1*Np1)/(d_i) + 4*Np1)*(rho_gm*u_g1^(02))/2;
d_Pt2s=(4*(1.58*log(Re_g2)-3.28)^(-02)*(L2*Np2)/(d_i) + 4*Np2)*(rho_gm*u_g2^(02))/2;
d_Ps1s=((exp(0.576-0.19*log(Re_s1)))*G_s1^(02)*(Nb1+1)*D1)/(2*rho_mlom*D_e1*phi_s1);
d_Ps2s=((exp(0.576-0.19*log(Re_s2)))*G_s2^(02)*(Nb2+1)*D2)/(2*rho_slom*D_e2*phi_s2);