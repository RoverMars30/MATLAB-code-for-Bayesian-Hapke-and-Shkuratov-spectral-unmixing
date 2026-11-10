function ssa_model = Hapke_forward_model(x, n_list, k_list, wl, isSpherical)
%--------------------------------------------------------------------------
% Hapke_forward_model:
%   根据输入端元参数计算混合 SSA（Hapke 模型）
%
% 输入：
%   x       - [2M x 1]，前 M 是 fraction，后 M 是粒径（μm）
%   n_list  - [L x M] 实部折射率
%   k_list  - [L x M] 虚部折射率
%   wl      - [L x 1] 波长（μm）
%   isSpherical - 是否球形颗粒
%
% 输出：
%   ssa_model - [L x 1] 混合 SSA
%--------------------------------------------------------------------------

    [L, M] = size(n_list);

    % -------- 拆参数 --------
    f = x(1:M);
    D = x(M+1:end);

    % f = f(:);
    % 
    % % -------- 可选：归一化（强烈推荐）--------
    % f(f < 0) = 0;
    % if sum(f) == 0
    %     f = ones(M,1)/M;
    % else
    %     f = f / sum(f);
    % end

    % -------- 初始化 --------
    ssa_model = zeros(L, 1);

    % -------- Hapke 混合 --------
    for j = 1:M
        n = n_list(:, j);
        k = k_list(:, j);

        ssa_j = Hapke_nk_to_ssa(wl, n, k, D(j), isSpherical);
        ssa_model = ssa_model + f(j) * ssa_j;
    end

end