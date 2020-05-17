function[xF,yF,psi]=flashisotermico(Tentrada,Pentrada,Tcentrada,Pcentrada,omegaentrada,ncentrada,Acentrada,fentrada)
%salidas=[xF,yF,psi]
%entradas=T,P,Tc,Pc,factor acentrico,nc
%constantes de antoine y flujos de entrada por componente
%T,P,Tc,Pc,omega,nc,Ac,f
%isotermico
nc=ncentrada;       %numero de componentes
Pc=Pcentrada;   %psi
Tc=Tcentrada;  %°F
omega=omegaentrada; 
Ac=Acentrada;
P=Pentrada;       %psi
f=fentrada;
F=sum(f);
z=f./F;

%Recibe T
Tvariable=Tentrada;

%supongo y o los calculo como ideales
ysup=[0.35;0.25;0.5;0.5;0.30];
toly=0.0001;
contadory=0;
cicloy=0;
while cicloy==0

%[keq]=calculodexyconTdefinido(Tvariable,P,Ac,Pc);
[keq]= SRK(ysup,P,Tvariable,nc,Pc,omega,Tc,Ac);
%yentrada,Pentrada,Tentrada,ncentrada,Pcentrada,omegaentrada,Tcentrada,Acentrada
psi=0.1;
tolerancia=0.001;
contador=0;
psical=0;
ep=1;

%utilizo Newton-Raphson para hallar psi 
while ep>tolerancia   
       
    fdepsi=sum((z.*(1-keq))./(1+psi.*(keq-1)));
    dfdepsi=sum((z.*(keq-1).^2)./(1+psi.*(keq-1)).^2);
    psical=psi-(fdepsi/dfdepsi);
    ep=abs((psical-psi)/psi);
    psi=psical;
    
    contador=contador+1;
    if contador>100
        contador
        'MUCHAS ITERACIONES'
        break
    end

end
psi
if psi>1
    psi=1;
elseif psi<0
    psi=0;
end

xFSN=z./(1+psi.*(keq-1));
xF=xFSN/sum(xFSN);
yFSN=keq.*xF;
yFcal=yFSN/sum(yFSN);
ep=abs((yFcal-ysup)./ysup);
epa=0;
for i=1:nc
    if ep(i)<toly
       epa=epa+1; 
    else
        epa=epa;
    end
end

if epa==nc
    cicloy=1;
    xF;
    yF=yFcal;
else
    ysup=yFcal;
   
end

contadory=contadory+1;
if contadory==50
    contadory
    'entró aqui'
    break
end

end

end
