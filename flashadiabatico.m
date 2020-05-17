%adiabatico
nc=3;       %numero de componentes
Pc=[714.2;587.8;530];   %psi
Tc=[1012.7;1069.1;1138];  %°F
omega=[0.2116;0.2415;0.2904]; 
format short e; bcp=[16.39282,4.02E-02,6.93E-06,-4.11E-08,2.40E-11;21.17722,4.64E-02,9.96E-06,-4.63E-08,2.59E-11;27.89247,5.10E-02,5.91E-06,-3.66E-08,1.95E-11];
format long; Ac=[5.658375,5307.813,379.456;5.944251,5836.287,374.745;5.922098,6141.641,354.0417];
format short;
P=14.7;       %psi
LK=1;
HK=3;       %Cuando son adyacentes
f=[136.077,181.436,136.077];
F=sum(f);
z=f./F;
Ru=1.987; %BTU/lbmol*R

%°F SUPUESTA
TF=164.0483;
Tb=(Ac(:,2)./(Ac(:,1)-log(P./Pc)))-Ac(:,3)
Tvector=[max(Tb),min(Tb),170];
%(max(Tb)+min(Tb))/2
for k=1:3
    
T=Tvector(k);     
[xF,yF,psi]=flashisotermico(T,P,Tc,Pc,omega,nc,Ac,f);
V=psi*F
L=F-V
Tref=TF; %Alimento como liquido saturado
integral=0;
for i=1:nc
    integral=integral+bcp(:,i).*((T^i-Tref^i)/i);
end
landa=(Ac(:,2).*Ru*(T+459.67)^2)./(Ac(:,3)+T).^2;
hV=integral+landa;
[correcHSRK]= SRKH(yF,P,T,nc,Pc,omega,Tc,Ac);
HV=sum(yF.*hV)+correcHSRK;
hL=integral;
HL=sum(xF.*hL);
HF=sum(z.*integral)
%HF=psi*HV+(1-psi)*HL
fdepsiconH=((psi*HV)+((1-psi)*HL)-HF)/1000;
fvector(k)=fdepsiconH;
end
Tvector
fvector




