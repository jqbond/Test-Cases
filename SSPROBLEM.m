function [] = l1
  k1f = 1000000;
  K1  = 30;
  k1r = k1f/K1;
  k2f = 1e-12;
  K2  = 1e-10;
  k2r = k2f/K2;
  k3f = 1e-12;
  K3  = 1;
  k3r = k3f/K3;

  param.k1f = k1f;
  param.k1r = k1r;
  param.k2f = k2f;
  param.k2r = k2r;
  param.k3f = k3f;
  param.k3r = k3r;
 
  ODEFUNC = @(t,var)(ODESYSTEM(t,var,param));
  var0 = [1;0;0;1];
  tspan = [0 1e5];
  [T,V] = ode15s(ODEFUNC, tspan, var0);
  thetaV = 1 - sum(V(:,2:3),2);

  figure(1)
  semilogx(T,V(:,2))
  legend('thetaA')

  figure(2)
  semilogx(T,V(:,3))
  legend('thetaB')
end

function [D] = ODESYSTEM(t,var, param)

  CA     = var(1);
  thetaA = var(2);
  thetaB = var(3);
  CB     = var(4);

  thetav = 1 - thetaA - thetaB;

  k1f = param.k1f;
  k1r = param.k1r;
  k2f = param.k2f;
  k2r = param.k2r;
  k3f = param.k3f;
  k3r = param.k3r;

  r1 = k1f*CA*thetav - k1r*thetaA;
  r2 = k2f*thetaA    - k2r*thetaB;
  r3 = k3f*thetaB    - k3r*CB;

  RA  =  0;
  RtA =  r1 - r2;
  RtB =  r2 - r3;
  RB  =  0;

  D = zeros(4,1);
  D(1) = RA;
  D(2) = RtA;
  D(3) = RtB;
  D(4) = RB;
#  [log10(abs(D(2:3))), log10(abs(var(2:3))), log10(abs(D(2:3).*var(2:3))), log10([t;t])]
#  pause(0.2) 
end
