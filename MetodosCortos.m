% Carga de bases de datos de componentes (se pueden anadir los que quiera,
% de manera consistente con las unidades)
load('bcp.mat');
load('Ac.mat');
load('omega.mat');
load('Pc.mat');
load('Tc.mat');
% Carga del formulario de command window
prompt = 'Escriba el numero de componentes: ';
nc = input(prompt);
prompt2 = 'Escriba la presion de la columna: ';
P = input(prompt2);
prompt2_1 = 'Escriba la temperatura de la alimentacion: ';
Tf = input(prompt2_1);
prompt2_2 = 'Escriba la presion de la alimentacion: ';
Pf = input(prompt2_2);
NoIdentificador = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'};
Componente = {'Etano';'Butano';'Acetona';'Hexano';'Acetato de etilo';'Benceno';'Heptano';'Octano';'o-Xileno';'Cumeno'};
Tt = table(NoIdentificador, Componente);
disp(Tt)
prompt3 = 'Escriba el flujo del primer componente: ';
ax = input(prompt3);
prompt4 = 'Escriba el flujo del segundo componente: ';
bx = input(prompt4);
prompt5 = 'Escriba el flujo del tercer componente: ';
cx = input(prompt5);
prompt6 = 'Escriba el flujo del cuarto componente: ';
dx = input(prompt6);
prompt7 = 'Escriba el flujo del quinto componente: ';
ex = input(prompt7);
prompt8 = 'Escriba el flujo del sexto componente: ';
fx = input(prompt8);
prompt9 = 'Escriba el flujo del septimo componente: ';
gx = input(prompt9);
prompt10 = 'Escriba el flujo del octavo componente: ';
hx = input(prompt10);
prompt11 = 'Escriba el flujo del noveno componente: ';
ix = input(prompt11);
prompt12 = 'Escriba el flujo del decimo componente: ';
jx = input(prompt12);
prompt13 = 'Escriba el numero identificador del clave ligero: ';
LK = input(prompt13);
prompt14 = 'Escriba la fraccion de recuperacion del clave ligero: ';
fLK = input(prompt14);
prompt15 = 'Escriba el numero identificador del clave pesado: ';
HK = input(prompt15);
prompt16 = 'Escriba la fraccion de recuperacion del clave pesado: ';
fHK = input(prompt16);
format short e; bcp;
format long; Ac;
format short;
f=[ax;bx;cx;dx;ex;fx;gx;hx;ix;jx];
F=sum(f);
z=f./F;
Ru=0.08314; %bar*m3/kmol*K
%estimado inical de las composiciones del destilado y fondos
d=[50;150;9.999133;99.37254;97;89.96542;1.5;0.00123345;0;0];
b=f-d;
tolciclogrande=0.001;
dios=0;
contadorciclogrande=0;

%Inicio de iteración
while dios==0
b=f-d;
F=sum(f); D=sum(d); B=sum(b);
xD=d./D;
xB=b./B;
z=f./F;

%****Calculo las temperaturas de ebullicion de los compuestos****
 Tb=Ac(:,2)./(Ac(:,1)-log(Pf./Pc))-Ac(:,3);
 Tmin=Tb(1);
 Tmax=Tb(nc);
 for i=1:nc
     if Tb(i)<Tmin
         Tmin=Tb(i);
     elseif Tb(i)>Tmax
         Tmax=Tb(i);
     else
     end
 end
 %****************fin temperaturas de ebullicion******************



%TEMPERATURA DEL DESTILADO************************************
Tvector=[Tmax,Tmin,(Tmax+Tmin)/2];
tolT=0.01;
seguridadT=0;
angel=0;
while angel==0
yD= z/3;
%yD=[0.45;0.25;0.3];  %Supuesto
for i=1:3
    Tvariable=Tvector(i);
    toly=0.001;
    ep=1;
    diablo=0;
    seguridad=0;
    while diablo==0
    
    [keq]=SRK(yD,P,Tvariable,nc,Pc,omega,Tc,Ac);
    yDcal=keq.*xD;
    yDN=yDcal./sum(yDcal);
    %ep=abs((yDN(LK)-yD(LK))/yD(LK));
    ep=abs((yDN-yD)./yD);
    epa=0;
    for j=1:nc
        if  ep(i)<toly
        epa=epa+1;
        end
    end
    if  epa==nc
        diablo=1;
        yD=yDN;
    else
      yD=yDN;
    end
    %seguridad
    seguridad=seguridad+1;
    if seguridad>50
        seguridad
        'Muchas iteraciones'
        break
    end
    end
    funvector(i)=sum(xD.*keq);
end


epT=abs(funvector(3)-1);
if epT<tolT
    TD=Tvector(3);
    angel=1;
elseif funvector(3)<1 && epT>tolT
    Tvector=[Tvector(1),Tvector(3),(Tvector(3)+Tvector(1))/2];
elseif funvector(3)>1 && epT>tolT
    Tvector=[Tvector(3),Tvector(2),(Tvector(3)+Tvector(2))/2];
end

seguridadT=seguridadT+1;
if seguridadT>50
        seguridadT
        'Muchas iteraciones de T'
        break
end
end
Keq(:,1)=keq;

%TEMPERATURA DEL FONDO***************************************
Tvector=[Tmax,Tmin,(Tmax+Tmin)/2];
tolT=0.01;
seguridadT=0;
angel=0;
while angel==0

yB= z/5;
%yB=[0.45;0.25;0.3;0;0];  %Supuesto
for i=1:3
    Tvariable=Tvector(i);
    toly=0.001;
    ep=1;
    diablo=0;
    seguridad=0;
    while diablo==0
    
    [keq]=SRK(yB,P,Tvariable,nc,Pc,omega,Tc,Ac);
    yBcal=keq.*xB;
    yBN=yBcal./sum(yBcal);
    %ep=abs((yBN(LK)-yB(LK))/yB(LK));
    ep=abs((yDN-yD)./yD);
    epa=0;
    for j=1:nc
        if  ep(i)<toly
        epa=epa+1;
        end
    end
    if  epa==nc
        diablo=1;
        yB=yBN;
    else
      yB=yBN;
    end
    %seguridad
    seguridad=seguridad+1;
    if seguridad>50
        seguridad
        'Muchas iteraciones'
        break
    end
    end
    
    funvector(i)=sum(xB.*keq);
end

epT=abs(funvector(3)-1);
if epT<tolT
    TB=Tvector(3);
    angel=1;
elseif funvector(3)<1 && epT>tolT
    Tvector=[Tvector(1),Tvector(3),(Tvector(3)+Tvector(1))/2];
elseif funvector(3)>1 && epT>tolT
    Tvector=[Tvector(3),Tvector(2),(Tvector(3)+Tvector(2))/2];
end

seguridadT=seguridadT+1;
if seguridadT>50
        seguridadT
        'Muchas iteraciones de T';
        break
end
end
Keq(:,2)=keq;

%TEMPERATURA DEL ALIMENTO***************************************
Tvector=[Tmax,Tmin,(Tmax+Tmin)/2];
tolT=0.01;
seguridadT=0;
angel=0;
while angel==0
yF = z;
%yF=[0.45;0.25;0.3;0;0];  %Supuesto
for i=1:3
    Tvariable=Tvector(i);
    toly=0.001;
    ep=1;
    diablo=0;
    seguridad=0;
    while diablo==0
    
    [keq]=SRK(yF,P,Tvariable,nc,Pc,omega,Tc,Ac);
    yFcal=keq.*z;
    yFN=yFcal./sum(yFcal);
    %ep=abs((yFN(LK)-yF(LK))/yF(LK));
    ep=abs((yDN-yD)./yD);
    epa=0;
    for j=1:nc
        if  ep(i)<toly
        epa=epa+1;
        end
    end
    if  epa==nc
        diablo=1;
        yF=yFN;
    else
      yF=yFN;
    end
    %seguridad
    seguridad=seguridad+1;
    if seguridad>50
        seguridad
        'Muchas iteraciones'
        break
    end
    end
    
    funvector(i)=sum(z.*keq);
end

epT=abs(funvector(3)-1);
if epT<tolT
    TF=Tvector(3);
    angel=1;
elseif funvector(3)<1 && epT>tolT
    Tvector=[Tvector(1),Tvector(3),(Tvector(3)+Tvector(1))/2];
elseif funvector(3)>1 && epT>tolT
    Tvector=[Tvector(3),Tvector(2),(Tvector(3)+Tvector(2))/2];
end

seguridadT=seguridadT+1;
if seguridadT>50
        seguridadT
        'Muchas iteraciones de T'
        break
end
end
Keq(:,3)=keq;
Keq; %Kde equilibrio de [Destilado, Fondo, Alimento]

TEMP=[TD,TB,TF];

for i=1:3
    ALPHA(:,i)=Keq(:,i)./Keq(HK,i);
end

alpham=(ALPHA(:,1).*ALPHA(:,2)).^0.5;
Nmin=log((d(LK)*b(HK))/(d(HK)*b(LK)))/log(alpham(LK));

bcal=f./(1+((d(HK)/b(HK))*(alpham).^Nmin));

dcal=f-bcal;

epciclogrande=abs((dcal-d)./d);
epa=0;

 for j=1:nc
        if  epciclogrande(j)<tolciclogrande
          epa=epa+1;
          epciclogrande(j);
        else
          epa=epa;  
        end
 end
 
    if  epa==nc
        epciclogrande;
        dios=1;
        d=dcal;
    else
        d=dcal;
    end
contadorciclogrande=contadorciclogrande+1;
if contadorciclogrande==100
    'Muchas iteraciones ciclo grande'
    break
end
end

'temperaturas de los flujos (TD,TB,TF)'
TEMP=[TD,TB,TF]
'Etapas minimas'
Nmin
'composiciones calculadas de los flujos(F  D  B)'
COMP=[f./F,d./D,b./B]
'Valores de K (KD,KB,KF)'
Keq
clase=(D.*xD)./(F.*z);
contadorclase=0;
for i=1:nc
if clase(i)>0 && clase(i)<1
    contadorclase=contadorclase+1;
else
    contadorclase=contadorclase;
end

end

%[xF,yF,psi]=flashisotermico(Tf,Pf,Tc,Pc,omega,nc,Ac,f)
%q=1-psi;
q=0.7;
LF=q*F;
if contadorclase~=nc
    'Clase 1'
    LinfsobreF=((LF/F)*((D*xD(LK)/(LF*z(LK)))-ALPHA(LK,3)*(D*xD(HK)/(LF*z(HK)))))/(ALPHA(LK,3)-1);
    Linf=LinfsobreF*F;
    %Para flujo molar constante
    'R minimo'
    Rminext=Linf/D
    
else 
    'Clase 2'
    %utilizo Newton-Raphson para hallar la raiz de Z del vapor 
 if HK-LK>=1 %**********CASO ADYACENTES*************
    theta=1.15;
    tolerancia=0.001;
    contador=0;
    thetacal=0;
    ep=1;

    while ep>tolerancia   
      
        fdetheta=sum((ALPHA(:,3).*z)./(ALPHA(:,3)-theta))-1+q;
        dfdetheta=sum((ALPHA(:,3).*z)./(ALPHA(:,3)-theta).^2);
        thetacal=theta-(fdetheta/dfdetheta);
        ep=abs((thetacal-theta)/theta);
        theta=thetacal;
    
        contador=contador+1;
        if contador>100
        contador
        'MUCHAS ITERACIONES'
        break
        end

    end
    
    'R minimo'
    Rminext=sum((ALPHA(:,3).*xD)./(ALPHA(:,3)-theta))-1
    %*************FIN ADYACENTES****************
 elseif HK-LK<1
        var=LK;
        for k=1:(HK-LK)
           theta=(ALPHA(var,3)+ALPHA(var+1,3))/2;
           %*********inicio Newton-Raphson***************
            tolerancia=0.001;
            contador=0;
            thetacal=0;
            ep=1;

            while ep>tolerancia   
    
            fdetheta=sum((ALPHA(:,3).*z)./(ALPHA(:,3)-theta))-1+q;
            dfdetheta=sum((ALPHA(:,3).*z)./(ALPHA(:,3)-theta).^2);
            thetacal=theta-(fdetheta/dfdetheta);
            ep=abs((thetacal-theta)/theta);
            theta=thetacal
    
            contador=contador+1;
            if contador>100
            contador
            'MUCHAS ITERACIONES'
            break
            end
            end
            
           %*********fin Newton-Raphson***************
           Raices(k)=theta;
           var=var+1;
           
        end
        
        
 end

    
end

'R operación'
Rop=1.2*Rminext
X=(Rop-Rminext)/(Rop+1);
Y=1-exp(((1+54.4*X)/(11+117.2*X))*((X-1)/(X^0.5)));
'Etapas reales'
Net=(Nmin+Y)/(1-Y)

%Calculo etapa de alimentación
NRsobreNS=((z(HK)/z(LK))*(xB(LK)/xD(HK))^2*(B/D))^0.206;
NS=(Net-1)/(NRsobreNS+1);
NR=NRsobreNS*NS;
'Etapa de alimentacion'
Alimentacion=NS+1

%Calculo de servicios

%encuentro la TEMPERATURA DE ROCIO del destilado 


Tvector=[Tmax,Tmin,(Tmax+Tmin)/2];
tolT=0.01;
seguridadT=0;
angel=0;
while angel==0
  yD=xD;
  for i=1:3
        Tvariable=Tvector(i);
        [keq]=SRK(yD,P,Tvariable,nc,Pc,omega,Tc,Ac);
        xDroc=yD./keq;
  end
    funvector(i)=sum(yD./keq);

    epT=abs(funvector(3)-1);
    if epT<tolT
         TDroc=Tvector(3);
         angel=1;
    elseif funvector(3)>1 && epT>tolT
         Tvector=[Tvector(1),Tvector(3),(Tvector(3)+Tvector(1))/2];
    elseif funvector(3)<1 && epT>tolT
        Tvector=[Tvector(3),Tvector(2),(Tvector(3)+Tvector(2))/2];
    end

seguridadT=seguridadT+1;
if seguridadT>50
        seguridadT
        'Muchas iteraciones de T'
        break
end
end


Tref=TDroc

landa=(Ac(:,2).*Ru*(TD+459.67)^2)./(Ac(:,3)+TD).^2;

integral=0;
for i=1:10
    integral = integral+((20*(TD+273.15) - 5463)*(5463*(TD+273.15)^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 20*(TD+273.15)^3*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 40*bcp(i,2).*bcp(i,3).^2*cosh(bcp(i,5)./(TD+273.15))^2 + 40*bcp(i,4).*bcp(i,5).^2*sinh(bcp(i,3)./(TD+273.15))^2))/(800*TD^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2);
end
integral;


hV=0;      %estado de referencia
hL=integral-landa;
HD=sum(hL.*xD)      %BTU/lbmol
'Calor del condensador'
QD=-(Rop+1)*D*HD     %BTU

%estado de referencia entalpia del destilado en estado liquido

HD=0;
Tref=TD;
integral1=0;

for i=1:10
    %integral1=integral1+bcp(:,i).*((TF^i-Tref^i)/i);
        integral1 = integral1+((20*(TD+273.15) - 5463)*(5463*(TD+273.15)^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 20*(TD+273.15)^3*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 40*bcp(i,2).*bcp(i,3).^2*cosh(bcp(i,5)./(TD+273.15))^2 + 40*bcp(i,4).*bcp(i,5).^2*sinh(bcp(i,3)./(TD+273.15))^2))/(800*TD^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2);

end
hF=integral1;
HF=sum(z.*hF);

Tref=TD;
integral2=0;
for i=1:10
    %integral2=integral2+bcp(:,i).*((TB^i-Tref^i)/i);
     integral2 = integral2+((20*(TD+273.15) - 5463)*(5463*(TD+273.15)^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 20*(TD+273.15)^3*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2 + 40*bcp(i,2).*bcp(i,3).^2*cosh(bcp(i,5)./(TD+273.15))^2 + 40*bcp(i,4).*bcp(i,5).^2*sinh(bcp(i,3)./(TD+273.15))^2))/(800*TD^2*cosh(bcp(i,5)./(TD+273.15))^2*sinh(bcp(i,3)./(TD+273.15))^2);

end
hB=integral2;
HB=sum(xB.*hB);
'Calor del ebullidor'
QB=D*HD+B*HB+QD-F*HF    %BTU

ALPHA
D
F
B

