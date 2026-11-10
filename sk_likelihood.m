function L = sk_likelihood(theta,D,LLK_prior)
    % theta: [2M x 1], 前 M 是粒径，后 M 是质量丰度
    % D: 数据结构体

    M = length(theta)/2;
    sizes = theta(1:M);
    wtperc = theta(M+1:end);

    fi = wtperctof(D.rho, sizes, wtperc);
    A_model = Shkuratov_forward_model([fi; sizes], D.n, D.k, D.lambda, D.q, D.isSpherical);

    res = D.R_SPEC - A_model;
    L = -0.5 * res' * D.Cinv * res;
end
