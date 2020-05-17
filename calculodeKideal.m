function [keq]= calculodeKideal(xentrada,Tentrada,Pentrada,Acentrada,Pcentrada)

Pc=Pcentrada;
Ac=Acentrada;
x=xentrada;

T=Tentrada;
psatantoine=Pc.*exp(Ac(:,1)-(Ac(:,2)./(Ac(:,3)+T)));
psat=psatantoine;
P=Pentrada;
alfa=0.1182+64.24/T;
beta=0.1735-42.27/T;
theta=0.1081;

gamma=exp(x.^2.*(alfa+beta.*(4*(1-x)-1)-theta*(1-2*x).*(5-6*x)));

keq=(gamma.*psat)./P;
 


end 