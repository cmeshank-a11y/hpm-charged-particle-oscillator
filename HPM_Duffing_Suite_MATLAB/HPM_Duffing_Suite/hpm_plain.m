function [X, dX, terms] = hpm_plain(eps_, tau, order)
%HPM_PLAIN  Plain HPM partial sum  X^(order) = sum_{n=0}^{order} eps^n u_n(tau).
%
%  [X, dX] = hpm_plain(eps_, tau, order) returns the partial sum and its
%  derivative dX/dtau.  terms(n+1,:) holds eps^n*u_n(tau) for n = 0..order.
%  Linear operator L = d2/dtau2 + 1, initial guess u_0 = cos(tau); this is
%  the formulation used in the original manuscript.  order <= 4.
if order > 4, error('hpm_plain: order must be <= 4'); end
U   = hpm_coefficients();
tau = tau(:).';
terms  = zeros(order+1, numel(tau));
dterms = zeros(order+1, numel(tau));
for r = 1:size(U,1)
    n = U(r,1);
    if n > order, continue; end
    c = U(r,2)/U(r,3); k = U(r,4); m = U(r,5); typ = U(r,6);
    tk  = tau.^k;
    if k > 0, dtk = k*tau.^(k-1); else, dtk = zeros(size(tau)); end
    if typ == 0
        f  = tk.*cos(m*tau);  df = dtk.*cos(m*tau) - m*tk.*sin(m*tau);
    else
        f  = tk.*sin(m*tau);  df = dtk.*sin(m*tau) + m*tk.*cos(m*tau);
    end
    terms(n+1,:)  = terms(n+1,:)  + eps_^n * c * f;
    dterms(n+1,:) = dterms(n+1,:) + eps_^n * c * df;
end
X  = sum(terms, 1);
dX = sum(dterms, 1);
end
