function [RP_result,GLL_result]=get_GLLintegral(D)

    lambda_all = D.lambda_all;
    mu_all = D.mu_all;
    rho_all = D.rho_all;
    kappa_all = D.kappa_all;

    xixl_all = D.xixl_all;
    xiyl_all = D.xiyl_all;
    xizl_all = D.xizl_all;
    etaxl_all = D.etaxl_all;
    etayl_all = D.etayl_all;
    etazl_all = D.etazl_all;

    gammaxl_all = D.gammaxl_all;
    gammayl_all = D.gammayl_all;
    gammazl_all = D.gammazl_all;

    wxgll_all = D.wxgll_all;
    wygll_all = D.wygll_all;
    wzgll_all = D.wzgll_all;

    jacob_all = D.jacob_all;

    points1 = D.points1;
    points2_reordered = D.points2_reordered;
    face_results = D.face_results;

    J_all_test = zeros(size(xixl_all));

    % Initialize result
    equal_indices = [];

    % Compare arrays element-wise to detect identical values
    for i = 1:length(xixl_all)
        if xixl_all(i) == xiyl_all(i) && xiyl_all(i) == xizl_all(i) && ...
           xizl_all(i) == etaxl_all(i) && etaxl_all(i) == etayl_all(i) && ...
           etayl_all(i) == etazl_all(i) && etazl_all(i) == gammaxl_all(i) && ...
           gammaxl_all(i) == gammayl_all(i) && gammayl_all(i) == gammazl_all(i)
            equal_indices = [equal_indices; i];
        end
    end

    % Zero-out duplicated derivative values
    xiyl_all(equal_indices)=0;
    xizl_all(equal_indices)=0;

    etaxl_all(equal_indices)=0;
    etazl_all(equal_indices)=0;

    gammaxl_all(equal_indices)=0;
    gammayl_all(equal_indices)=0;

    % Compute inverse Jacobian determinant at each point
    for i = 1:length(xixl_all)

        xixl = xixl_all(i);
        xiyl = xiyl_all(i);
        xizl = xizl_all(i);

        etaxl = etaxl_all(i);
        etayl = etayl_all(i);
        etazl = etazl_all(i);

        gammaxl = gammaxl_all(i);
        gammayl = gammayl_all(i);
        gammazl = gammazl_all(i);

        J_inv = calculate_inverse_jacobian(xixl, xiyl, xizl, etaxl, etayl, etazl, ...
                                           gammaxl, gammayl, gammazl, 1);

        J_all_test(i) = det(J_inv);
    end

    U_RP=0;

    % ------------------------------------------------------
    % Compute surface integrals over each boundary
    % ------------------------------------------------------

    % ------------------- pxmin surface --------------------
    U_RP_all = face_results.pxmin.U_RP_all;
    U_RP_pxmin = sum(U_RP_all,2);
    pxmin = face_results.pxmin.faceidx;

    [weights_matrix_pxmin,pxmin_integral] = compute_surface_integral(...
         U_RP_all, xixl_all, xiyl_all, xizl_all, ...
         etaxl_all, etayl_all, etazl_all, ...
         gammaxl_all, gammayl_all, gammazl_all, ...
         wxgll_all, wygll_all, wzgll_all, pxmin, 'yz');

    U_RP = U_RP + U_RP_pxmin;

    % ------------------- pxmax surface --------------------
    U_RP_all = face_results.pxmax.U_RP_all;
    U_RP_pxmax = sum(U_RP_all,2);
    pxmax = face_results.pxmax.faceidx;

    [weights_matrix_pxmax,pxmax_integral] = compute_surface_integral(...
         U_RP_all, xixl_all, xiyl_all, xizl_all, ...
         etaxl_all, etayl_all, etazl_all, ...
         gammaxl_all, gammayl_all, gammazl_all, ...
         wxgll_all, wygll_all, wzgll_all, pxmax, 'yz');

    U_RP = U_RP + U_RP_pxmax;

    % ------------------- pymin surface --------------------
    U_RP_all = face_results.pymin.U_RP_all;
    U_RP_pymin = sum(U_RP_all,2);
    pymin = face_results.pymin.faceidx;

    [weights_matrix_pymin,pymin_integral] = compute_surface_integral(...
         U_RP_all, xixl_all, xiyl_all, xizl_all, ...
         etaxl_all, etayl_all, etazl_all, ...
         gammaxl_all, gammayl_all, gammazl_all, ...
         wxgll_all, wygll_all, wzgll_all, pymin, 'xz');

    U_RP = U_RP + U_RP_pymin;

    % ------------------- pymax surface --------------------
    U_RP_all = face_results.pymax.U_RP_all;
    U_RP_pymax = sum(U_RP_all,2);
    pymax = face_results.pymax.faceidx;

    [weights_matrix_pymax,pymax_integral] = compute_surface_integral(...
         U_RP_all, xixl_all, xiyl_all, xizl_all, ...
         etaxl_all, etayl_all, etazl_all, ...
         gammaxl_all, gammayl_all, gammazl_all, ...
         wxgll_all, wygll_all, wzgll_all, pymax, 'xz');

    U_RP = U_RP + U_RP_pymax;

    % ------------------- pzbottom surface -----------------
    U_RP_all = face_results.pzbottom.U_RP_all;
    U_RP_pzbottom = sum(U_RP_all,2);
    pzbottom = face_results.pzbottom.faceidx;

    [weights_matrix_pzbottom,pzbottam_integral] = compute_surface_integral(...
         U_RP_all, xixl_all, xiyl_all, xizl_all, ...
         etaxl_all, etayl_all, etazl_all, ...
         gammaxl_all, gammayl_all, gammazl_all, ...
         wxgll_all, wygll_all, wzgll_all, pzbottom, 'xy');

    U_RP = U_RP + U_RP_pzbottom;

    % Total GLL integral
    total_integral = pxmin_integral + pxmax_integral + ...
                     pymin_integral + pymax_integral + ...
                     pzbottam_integral;

    total_integral_sum = U_RP;

    % Store GLL results
    GLL_result.total_integral      = total_integral;
    GLL_result.pxmin_integral      = pxmin_integral;
    GLL_result.pxmax_integral      = pxmax_integral;
    GLL_result.pymin_integral      = pymin_integral;
    GLL_result.pymax_integral      = pymax_integral;
    GLL_result.pzbottam_integral   = pzbottam_integral;

    % Store RP results
    RP_result.total_integral_sum = total_integral_sum;
    RP_result.U_RP_pxmin = U_RP_pxmin;
    RP_result.U_RP_pxmax = U_RP_pxmax;
    RP_result.U_RP_pymin = U_RP_pymin;
    RP_result.U_RP_pymax = U_RP_pymax;
    RP_result.U_RP_pzbottom = U_RP_pzbottom;

% ----------------------------------------------------------
% Helper functions
% ----------------------------------------------------------

function J_inv = calculate_inverse_jacobian(xixl, xiyl, xizl, ...
         etaxl, etayl, etazl, gammaxl, gammy, gammz, index)
    % Construct 3×3 Jacobian matrix
    J = [xixl(index),    xiyl(index),    xizl(index);
         etaxl(index),   etayl(index),   etazl(index);
         gammaxl(index), gammy(index),   gammz(index)];

    % Compute inverse
    J_inv = inv(J);
end

function J_inv = calculate_2d_jacobian_inverse(xixl, xiyl, etaxl, etayl, index)
    % Construct 2D Jacobian matrix for XY plane
    J = [xixl(index), xiyl(index);
         etaxl(index), etayl(index)];

    J_inv = inv(J);
end

function J_inv = calculate_2d_jacobian_inverse_xz(xixl, xizl, ...
         gammaxl, gammazl, index)
    % Construct 2D Jacobian matrix for XZ plane
    J = [xixl(index), xizl(index);
         gammaxl(index), gammazl(index)];

    J_inv = inv(J);
end

function J_inv = calculate_2d_jacobian_inverse_yz(etayl, etazl, ...
         gammay, gammaz, index)
    % Construct 2D Jacobian matrix for YZ plane
    J = [etayl(index), etazl(index);
         gammay(index), gammaz(index)];

    J_inv = inv(J);
end

% ----------------------------------------------------------
% Compute surface GLL integral
% ----------------------------------------------------------
function [weights_matrix,area_integral] = compute_surface_integral(...
    Ux_RP_all, xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, point_indices, direction)

    [num_time_steps, num_space_points] = size(Ux_RP_all);

    % Initialize outputs
    area_integral = zeros(num_time_steps, 1);
    weights_matrix = zeros(length(point_indices), 1);

    for p = 1:length(point_indices)
        index = point_indices(p);

        % Select Jacobian type and GLL weights
        if strcmp(direction, 'yz')
            J_inv = calculate_2d_jacobian_inverse_yz(etayl_all, etazl_all, ...
                                                     gammayl_all, gammazl_all, index);
            weight = wygll_all(index) * wzgll_all(index);

        elseif strcmp(direction, 'xz')
            J_inv = calculate_2d_jacobian_inverse_xz(xixl_all, xizl_all, ...
                                                     gammaxl_all, gammazl_all, index);
            weight = wxgll_all(index) * wzgll_all(index);

        elseif strcmp(direction, 'xy')
            J_inv = calculate_2d_jacobian_inverse(xixl_all, xiyl_all, ...
                                                  etaxl_all, etayl_all, index);
            weight = wxgll_all(index) * wygll_all(index);
        end

        % Avoid NaN determinant
        if isnan(det(J_inv))
            J_det = 1 / 0.002^2;
        else
            J_det = det(J_inv);
        end

        area_integral = area_integral + Ux_RP_all(:, p) * J_det * weight;
        weights_matrix(p) = J_det * weight;
    end

end

end
