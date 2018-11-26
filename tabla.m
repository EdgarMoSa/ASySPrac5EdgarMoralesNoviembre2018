Dn = {'D0';'D1';'D2';'D3';'D4'};
Real = [0.504;0.0296 - 0.1186i;0.0078 - 0.0620i;0.0035 - 0.0417i;0.0020 - 0.0314i];
Trapecio = [0.5047+0.0000i;0.0301-0.1168i;0.0082-0.0583i;0.0040-0.0361i;0.0025-0.0237i];
Lathi = [0.504281105948516 + 0.00000000000000i;0.0296650836043681 - 0.118647676913452i;0.00775972886772522 - 0.0620525144086751i;0.00347937240545304 - 0.0417144902523740i;0.00196376005535414 - 0.0313695156147827i];
TrapeciovsReal = abs(Real-Trapecio);
LathivsReal = abs(Real-Lathi);
T=table(Dn,Real,Trapecio,Lathi,TrapeciovsReal,LathivsReal)