function [rb, rf] = Shkuratov_Calc_rb_rf(n, k, D, lambda, isSpherical)
    % n, k: (Lx1)
    % 返回每个波长的 rb, rf

    if isSpherical == 0
        S = 0.2 * D;
    else
        S = (2/3) * (n.^2 - (1./n) .* (n.^2 - 1).^(3/2)) .* D;
    end

    r0 = ((n - 1) ./ (n + 1)).^2;
    Re = r0 + 0.05;
    Rb = (0.28 .* n - 0.20) .* Re;
    Ri = 1.04 - 1 ./ (n.^2);
    Te = 1 - Re;
    Ti = Te ./ n.^2;

    tau = (4 * pi .* k .* S) ./ lambda;
    exp_tau = exp(-tau);
    exp_2tau = exp(-2 * tau);

    numerator = 0.5 * Te .* Ti .* Ri .* exp_2tau;
    denominator = 1 - Ri .* exp_tau;

    rb = Rb + numerator ./ denominator;
    rf = (Re - Rb) + Te .* Ti .* exp_tau + numerator ./ denominator;
end
