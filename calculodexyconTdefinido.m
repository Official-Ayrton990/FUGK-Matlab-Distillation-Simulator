function [keq]= calculodexyconTdefinido(Tentrada,Pentrada,Acentrada,Pcentrada)

Pc=Pcentrada;
Ac=Acentrada;
%z=zentrada;
P=Pentrada;
T=Tentrada;
psatantoine=Pc.*exp(Ac(:,1)-(Ac(:,2)./(Ac(:,3)+T)));
psat=psatantoine;
keq=psat./P;

end 