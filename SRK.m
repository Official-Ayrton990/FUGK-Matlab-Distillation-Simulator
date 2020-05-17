

%*****************SRK*****************
function [keq]=SRK(yentrada,Pentrada,Tentrada,ncentrada,Pcentrada,omegaentrada,Tcentrada,Acentrada)
%RECIBE: composicion del vapor,P(psi),T(°F),nc,Pc,factor acentrico,Tc y constantes de antoine

%CALCULA la k de equilibrio de los componentes
  

PRESION=Pentrada;   %psia
Pc=Pcentrada;   %psia
Tc=Tcentrada; %°R
factascen=omegaentrada;
nc=ncentrada;           %numero de componentes
Ac=Acentrada;
T=Tentrada;


y= yentrada;   %datos que ingresan por ahora los fijo


Psat=Pc.*exp(Ac(:,1)-(Ac(:,2)./(Ac(:,3)+T)));

% %Calculo las temperaturas de ebullicion
% pantoine=PRESION;     %Presion en psi
% Tbantoine=Ac(:,2)./(Ac(:,1)-log(pantoine./Pc))-Ac(:,3);
% Tb=(5/9).*(Tbantoine-32);
% %fin del calculo de T de ebullicion de las especies puras

Tr=((T+459.67))./Tc;
Pr=(PRESION)./Pc;

msrk=0.480+1.574.*factascen+0.176.*factascen.^2;
amin=(1+msrk.*(1-Tr.^0.5)).^2;

Aast=0.42747.*amin.*(Pr./Tr.^2);
Bast=0.08664.*(Pr./Tr);

Bmez=0;
Amez=0;
for i=1:nc
    Bmez=Bmez+y(i)*Bast(i);
    for j=1:nc
        Amez=Amez+y(i)*y(j)*(Aast(i)*Aast(j))^0.5;
    end
end

Z1=1;
tolerancia=0.001;
contador=0;
Z1cal=0;
ep=1;

%utilizo Newton-Raphson para hallar la raiz de Z del vapor 
while ep>tolerancia   
    
    
    fdeZ=Z1^3-Z1^2+Z1*(Amez-Bmez-Bmez^2)-Amez*Bmez;
    dfdeZ=3*Z1^2-2*Z1+Amez-Bmez-Bmez^2;
    Z1cal=Z1-(fdeZ/dfdeZ);
    ep=abs((Z1cal-Z1)/Z1);
    
    Z1=Z1cal;
    
    contador=contador+1;
    if contador>100
        contador
        'MUCHAS ITERACIONES'
        break
    end

end

Z=Z1;

fi=exp((Z-1)*(Bast./Bmez)-log(Z-Bmez)-(Amez/Bmez).*(2.*((Aast.^0.5)./Amez^0.5)-(Bast./Bmez)).*log((Z+Bmez)./Z));
keq=Psat./(fi*PRESION);

end
