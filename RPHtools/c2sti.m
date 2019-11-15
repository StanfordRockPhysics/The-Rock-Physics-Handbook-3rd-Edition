function [bb]=c2sti(aa)
% [BB]=C2STI(AA)
% convert between stiffness and compliance for TI (or hexagonal) symmetry.
% If the input is the stiffness, output is compliance and vice-versa.
%
% input:  AA=(a11,a12,a13,a33,a44)
% output: BB=(b11,b12,b13,b33,b44)

% Written by Frank Liu

a11=aa(:,1);
a12=aa(:,2);
a13=aa(:,3);
a33=aa(:,4);
a44=aa(:,5);

as=a33.*(a11+a12)-2*a13.*a13;

b11=(a33./as+1./(a11-a12))/2;
b12=(a33./as-1./(a11-a12))/2;
b13=-a13./as;
b33=(a11+a12)./as;
b44=1./a44;

bb=[b11 b12 b13 b33 b44];

