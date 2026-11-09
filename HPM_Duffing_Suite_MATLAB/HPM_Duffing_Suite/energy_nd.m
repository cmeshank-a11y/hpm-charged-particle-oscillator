function E = energy_nd(eps_, X, dX)
%ENERGY_ND  Non-dimensional Hamiltonian H = X'^2/2 + X^2/2 + eps*X^4/4.
%  The physical energy per unit mass is C*A^2*H, so relative variations are identical.
E = 0.5*dX.^2 + 0.5*X.^2 + 0.25*eps_*X.^4;
end
