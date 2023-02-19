function [A,B,C,Cz,E] = linearize(p, hs)
    ap = p(1:4);
    At = p(5:8);
    gam = p(9:10);
    g = p(11);
    rho = p(12);

    T = (At./ap).*sqrt(2*hs/g);

    A=[-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
    B=[rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
    C=diag(1./(rho*At));
    Cz=C(1:2,:);
    E = [0 0; 0 0; rho 0; 0 rho];
end

