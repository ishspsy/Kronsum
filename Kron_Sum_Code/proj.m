function[proj_inv]=proj(invS);

[si,si1]=size(invS);

cvx_begin quiet;
   variable x(si,si) symmetric;
   minimize(norm(x-invS,1));
cvx_end;

proj_inv=x;
