function [sigma] = sigParra(f1,f2,gama,zeta)

%sigma=1./(gama.^2).*(1-(sin(pi*f1*zeta).^2).*(sin(pi*f2*zeta).^2)./(pi^4*zeta.^4*f1^2*f2^2));
sigma=(sin(pi*f1*zeta/gama)).*(sin(pi*f2*zeta/gama))./(pi^2*zeta.^2*f1*f2);