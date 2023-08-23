function dydt = ddt_EVO(t,y)

global u ks k d eff A1 A2 Nm K

y(y<0) = 0;

dydt = zeros(length(y),1);
tmp = 0;
for ii = 1:length(y)-1
    % MAIN model
    % dNdt
    dydt(ii) = (1 - sum(y(1:end-1)) ./ Nm) .* u(ii).*y(ii) .* (y(end)./(y(end)+ks(ii))) .* k(ii)./(k(ii)+A1) - d(ii)*A2*y(ii)*y(end)./(ks(ii)+y(end))*eff(ii)*3e6;

    % dSdt
    tmp = tmp - eff(ii).*y(ii).*u(ii).*y(end)./(y(end)+ks(ii));
end

dydt(end) = tmp;

