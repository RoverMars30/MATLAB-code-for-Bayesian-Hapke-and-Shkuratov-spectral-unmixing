function A_model = Shkuratov_forward_model(x, n_list, k_list, lambda, q, isSpherical)
    % 根据输入参数 x（比例和S），计算模拟的反照率 A_model

    % 计算混合颗粒物反射率的 Shkuratov 模型（粗颗粒近似）
    % 输入：
    %   x       - [2M x 1]，前 M 是丰度，后 M 是粒径（μm）
    %   n_list  - [L x M] 折射率实部
    %   k_list  - [L x M] 折射率虚部
    %   lambda  - [L x 1] 波长（μm）
    %   q       - 填充率（0~1）
    % 输出：
    %   A_model - [L x 1] 模拟反射率

    L = length(lambda);
    total_len = length(x);
    M = total_len / 2;

    f = x(1:M);                 % 各端元比例
    S_vec = x(M+1:end);         % 每个端元的 S 值

    % 检查维度一致性
    if size(n_list, 1) ~= L || size(k_list, 1) ~= L
        error('n_list 和 k_list 的行数应等于波段数 L');
    end
    if size(n_list, 2) ~= M || size(k_list, 2) ~= M
        error('n_list 和 k_list 的列数应等于 M');
    end

    rb_all = zeros(L, M);
    rf_all = zeros(L, M);
    
    f = f(:);
    
    for j = 1:M
        n = n_list(:, j);
        k = k_list(:, j);
        S_j = S_vec(j);
        [rb, rf] = Shkuratov_Calc_rb_rf(n, k, S_j, lambda, isSpherical); % isSpherical = 1
        rb_all(:, j) = rb;
        rf_all(:, j) = rf;
    end

    rho_b = q * (rb_all * f);
    rho_f = q * (rf_all * f) + (1 - q);
    term = (1 + rho_b.^2 - rho_f.^2) ./ (2 * rho_b);
    A_model = term - sqrt(term.^2 - 1);
end