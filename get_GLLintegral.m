function [RP_result, GLL_result] = get_GLLintegral(D)
% Author: Jingnan Sun
% Affiliation: Peking University
% Date: 2025.11
% ========================================================================
% get_GLLintegral
%
% This function performs the Gauss–Lobatto–Legendre (GLL) surface
% integration for the Representation Theorem (RP) boundary terms produced
% by prepare_data.m. It:
%
%   1. Loads geometric operators and material parameters from D
%   2. Calculate Jacobian matrix
%   3. Integrates displacement (U_RP) contributions over all five boundary surfaces:
%          pxmin, pxmax, pymin, pymax, pzbottom
        %           pymax
        % ---------------------------------
%   pxmin |                               | pxmax
        % |                               |
        % |           pzbottom            |    z-axis ↑ (out of plane)
        % |                               |
        % |                               |
        % ---------------------------------
        %           pymin
        %            y ↑
        %           x →


%   4. Produces:
%         - GLL_result: final boundary integrals for each surface
%         - RP_result : raw RP term sums for backup check
%
% The output waveform is:
%       GLL_result.total_integral
%
% ========================================================================


% ------------------------------------------------------------------------
% Extract material parameters
% These material parameters are extracted from SPECFEM3D_Cartesian (Komatitsch & Tromp, 1999),
% and can be reconstructed from the files provided in the binfiles directory.
% ------------------------------------------------------------------------

lambda_all = D.lambda_all;
mu_all     = D.mu_all;
rho_all    = D.rho_all;
kappa_all  = D.kappa_all;



% ------------------------------------------------------------------------
% Extract derivatives in the natural coordinate system of spectral element.
% These define ∂ξ/∂x, ∂ξ/∂y, ∂ξ/∂z etc. for each GLL node
%       xixl  = ∂ξ/∂x
%       xiyl  = ∂ξ/∂y
%       xizl  = ∂ξ/∂z
%       etaxl = ∂η/∂x 
%       ...
%  coordinate transformation:
% 
% (1) From natural coordinates (ξ, η, γ) to physical coordinates (x, y, z):
%     [ x ]     [ ∂x/∂ξ  ∂x/∂η  ∂x/∂γ ]   [ ξ ]
%     [ y ]  =  [ ∂y/∂ξ  ∂y/∂η  ∂y/∂γ ] * [ η ]
%     [ z ]     [ ∂z/∂ξ  ∂z/∂η  ∂z/∂γ ]   [ γ ]
% 
%     Compact form:
%         x_phys = J * ξ_nat
% 
% (2) From physical coordinates (x, y, z) back to natural coordinates (ξ, η, γ):
%     [ ξ ]     [ ∂ξ/∂x  ∂ξ/∂y  ∂ξ/∂z ]   [ x ]
%     [ η ]  =  [ ∂η/∂x  ∂η/∂y  ∂η/∂z ] * [ y ]
%     [ γ ]     [ ∂γ/∂x  ∂γ/∂y  ∂γ/∂z ]   [ z ]
% 
%     Compact form:
%         ξ_nat = J^{-1} * x_phys
% 
% Relationship:
%         J^{-1} = ∂(ξ,η,γ)/∂(x,y,z)
%         J     = ∂(x,y,z)/∂(ξ,η,γ)

% see Komatitsch & Tromp (1999) for their detailed definitions.
% ------------------------------------------------------------------------
xixl_all = D.xixl_all;
xiyl_all = D.xiyl_all;
xizl_all = D.xizl_all;

etaxl_all = D.etaxl_all;
etayl_all = D.etayl_all;
etazl_all = D.etazl_all;

gammaxl_all = D.gammaxl_all;
gammayl_all = D.gammayl_all;
gammazl_all = D.gammazl_all;



% ------------------------------------------------------------------------
% Extract GLL weights and Jacobian
%
% 2-D Gauss–Lobatto–Legendre (GLL) quadrature rule:
%
%     ∫_{-1}^{1} ∫_{-1}^{1} f(ξ,η) dA  ≈  Σ_i Σ_j  w_i w_j f(ξ_i, η_j) * |J_ij|
%
% where:
%   ξ_i, η_j   = GLL nodes in [-1, 1]
%   w_i, w_j   = associated GLL quadrature weights
%   J_ij       = 2-D surface Jacobian determinant at node (i,j)
%
% For a 4th-order GLL quadrature (5 GLL points), the nodes are:
%
%     ξ = (-1,  -sqrt(3/7),  0,  sqrt(3/7),  1)
%
% and the weights are:
%
%     w = ( 9/90,  49/90,  64/90,  49/90,  9/90 )
%
% The grid of nodes and weights looks like this:
%
%
%          η = 1        49/90  •────────•────────•────────•────────• 49/90
%                              |        |        |        |       |
%          η = sqrt(3/7) 64/90 •────────•────────•────────•────────• 64/90
%                              |        |        |        |       |
%          η = 0        49/90  •────────•────────•────────•────────• 49/90
%                              |        |        |        |       |
%          η = -sqrt(3/7)64/90 •────────•────────•────────•────────• 64/90
%                              |        |        |        |       |
%          η = -1        49/90 •────────•────────•────────•────────• 49/90
%
%                        ξ = -1      -√(3/7)       0      √(3/7)      1
%
%
% At each node (ξ_i,η_j), the surface integral uses:
%
%       w_i * w_j * |J(ξ_i,η_j)|
% to represent the area of subsurface dA.
% SPECFEM3D_Cartesian provides the Jacobian components in inverse form
% (∂ξ/∂x, ∂η/∂x, ...), so the standard surface Jacobian must be obtained by
% inverting this matrix.
%
% see Komatitsch & Tromp (1999) for their detailed definitions.
% ------------------------------------------------------------------------
wxgll_all = D.wxgll_all;
wygll_all = D.wygll_all;
wzgll_all = D.wzgll_all;

jacob_all = D.jacob_all;   % (3-D Jacobian, not used for RP 2-D surfaces)


% ------------------------------------------------------------------------
% Load node coordinates and face results produced in prepare_data.m
% ------------------------------------------------------------------------
points1           = D.points1;
points2_reordered = D.points2_reordered;
face_results      = D.face_results;


% ------------------------------------------------------------------------
% Preallocate Jacobian determinant test array (For verification purposes only)
% ------------------------------------------------------------------------
J_all_test = zeros(size(xixl_all));


% ========================================================================
% STEP 1 — Detect repeated derivative fields and correct them
% Some derivative fields may be repeated due to regular elements.
% Those values are zeroed out to avoid producing singular Jacobians.
% ------------------------------------------------------------------------
% NOTE:
% This is not hard-coding. The behavior is required to match the file
% format generated by SPECFEM3D_Cartesian.
%
% In SPECFEM3D_Cartesian, to reduce storage size, regular hexahedral
% elements store only one derivative value per GLL point, because the
% three spatial derivatives (∂/∂ξ, ∂/∂η, ∂/∂γ) are identical at that node.
% in the regular element
%
% Therefore, repeated derivative fields appearing in the input data are a
% feature of the SPECFEM storage scheme, not a coding issue here.
% ------------------------------------------------------------------------
% ========================================================================

equal_indices = [];

for i = 1:length(xixl_all)
    if xixl_all(i) == xiyl_all(i) && ...
       xiyl_all(i) == xizl_all(i) && ...
       xizl_all(i) == etaxl_all(i) && ...
       etaxl_all(i) == etayl_all(i) && ...
       etayl_all(i) == etazl_all(i) && ...
       etazl_all(i) == gammaxl_all(i) && ...
       gammaxl_all(i) == gammayl_all(i) && ...
       gammayl_all(i) == gammazl_all(i)

        equal_indices = [equal_indices; i];
    end
end

% ------------------------------------------------------------------------
% Zero-out derivatives for regular elements
%
% NOTE:
% This is not hard-coding. This correction is required to match the way
% SPECFEM3D_Cartesian stores geometric derivatives for *regular elements*.
%
% In SPECFEM3D_Cartesian, for regular hexahedral elements the derivative
% fields are stored in a compressed form:
%
%       ξ_x = (∂ξ/∂x)_regular,   but  ξ_y and ξ_z are NOT stored separately
%
% Because for regular elements:
%
%       ∂ξ/∂x = ξ_x_regular
%       ∂ξ/∂y = 0
%       ∂ξ/∂z = 0
%
% and similarly:
%
%       ∂η/∂x = 0,   ∂η/∂y = η_y_regular,   ∂η/∂z = 0
%       ∂γ/∂x = 0,   ∂γ/∂y = 0,             ∂γ/∂z = γ_z_regular
%
% Because the SPECFEM output stores only the nonzero component, the other
% two components must be explicitly set to zero. These assignments restore
% the correct full derivative structure required for constructing the
% geometric Jacobian.
%

xiyl_all(equal_indices) = 0;
xizl_all(equal_indices) = 0;

etaxl_all(equal_indices) = 0;
etazl_all(equal_indices) = 0;

gammaxl_all(equal_indices) = 0;
gammayl_all(equal_indices) = 0;


% ========================================================================
% STEP 2 — Compute inverse Jacobian determinant for all nodes (For verification purposes only)
% Users can compare J_all_test (calculated by formula) 
% 3-D Jacobian matrix for SEM coordinate transformation:
%
%       [ ∂x/∂ξ   ∂x/∂η   ∂x/∂γ ]
%   J = [ ∂y/∂ξ   ∂y/∂η   ∂y/∂γ ]
%       [ ∂z/∂ξ   ∂z/∂η   ∂z/∂γ ]
%
% The Jacobian determinant is:
%
%   detJ = det(J)
%
% This transforms derivatives between natural coordinates (ξ,η,γ)
% and physical coordinates (x,y,z).

% with jacob_all (get directly from SPECFEM3D_Cartesian) to verify 

% ========================================================================
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

    J_inv = calculate_inverse_jacobian( ...
                xixl, xiyl, xizl, ...
                etaxl, etayl, etazl, ...
                gammaxl, gammayl, gammazl, 1 );

    J_all_test(i) = det(J_inv);
end


% ========================================================================
% STEP 3 — Initialize RP accumulator
% ========================================================================
U_RP = 0;


% ========================================================================
% STEP 4 — Surface-by-surface GLL integration
% Each surface contributes: pxmin, pxmax, pymin, pymax, pzbottom
% ========================================================================


% ------------------------------------------------------------------------
% PXMIN
% ------------------------------------------------------------------------
U_RP_all   = face_results.pxmin.U_RP_all;
U_RP_pxmin = sum(U_RP_all, 2);
pxmin      = face_results.pxmin.faceidx;

[weights_matrix_pxmin, pxmin_integral] = compute_surface_integral( ...
    U_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    pxmin, 'yz');

U_RP = U_RP + U_RP_pxmin;


% ------------------------------------------------------------------------
% PXMAX
% ------------------------------------------------------------------------
U_RP_all   = face_results.pxmax.U_RP_all;
U_RP_pxmax = sum(U_RP_all,2);
pxmax      = face_results.pxmax.faceidx;

[weights_matrix_pxmax, pxmax_integral] = compute_surface_integral( ...
    U_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    pxmax, 'yz');

U_RP = U_RP + U_RP_pxmax;


% ------------------------------------------------------------------------
% PYMIN
% ------------------------------------------------------------------------
U_RP_all   = face_results.pymin.U_RP_all;
U_RP_pymin = sum(U_RP_all,2);
pymin      = face_results.pymin.faceidx;

[weights_matrix_pymin, pymin_integral] = compute_surface_integral( ...
    U_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    pymin, 'xz');

U_RP = U_RP + U_RP_pymin;


% ------------------------------------------------------------------------
% PYMAX
% ------------------------------------------------------------------------
U_RP_all   = face_results.pymax.U_RP_all;
U_RP_pymax = sum(U_RP_all,2);
pymax      = face_results.pymax.faceidx;

[weights_matrix_pymax, pymax_integral] = compute_surface_integral( ...
    U_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    pymax, 'xz');

U_RP = U_RP + U_RP_pymax;


% ------------------------------------------------------------------------
% PZBOTTOM
% ------------------------------------------------------------------------
U_RP_all      = face_results.pzbottom.U_RP_all;
U_RP_pzbottom = sum(U_RP_all,2);
pzbottom      = face_results.pzbottom.faceidx;

[weights_matrix_pzbottom, pzbottam_integral] = compute_surface_integral( ...
    U_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    pzbottom, 'xy');

U_RP = U_RP + U_RP_pzbottom;


% ========================================================================
% STEP 5 — Combine total boundary integral
% ========================================================================
total_integral = pxmin_integral + pxmax_integral + ...
                 pymin_integral + pymax_integral + ...
                 pzbottam_integral;

total_integral_sum = U_RP;


% ========================================================================
% STEP 6 — Save GLL integral results
% ------------------------------------------------------------------------
% NOTE ABOUT THE DIFFERENCE BETWEEN GLL_result AND RP_result:
%
% GLL_result:
%   - These waveforms are computed using the *correct* GLL surface
%     integration:
%         • GLL quadrature weights (w_i, w_j)
%         • 2-D Jacobian determinant |J|
%         • Proper summation over all boundary nodes
%   - Therefore, GLL_result contains the physically meaningful,
%     fully-integrated boundary contributions for each surface.
%
% RP_result:
%   - These are diagnostic-only waveforms obtained by simply summing the
%     raw U_RP contributions over all nodes on each face *without* applying
%     any weights or Jacobian scaling.
%   - They represent an "equal-weight" sum of U_RP(t) at each node and are
%     NOT a proper numerical integral.
%
% In short:
%   • Use GLL_result for the correct physical synthetic waveform.
%   • RP_result is only for debugging.
% ------------------------------------------------------------------------
% ========================================================================
GLL_result.total_integral    = total_integral;
GLL_result.pxmin_integral    = pxmin_integral;
GLL_result.pxmax_integral    = pxmax_integral;
GLL_result.pymin_integral    = pymin_integral;
GLL_result.pymax_integral    = pymax_integral;
GLL_result.pzbottam_integral = pzbottam_integral;


% ========================================================================
% STEP 7 — Save RP results (just for backup)
% ========================================================================
RP_result.total_integral_sum = total_integral_sum;
RP_result.U_RP_pxmin         = U_RP_pxmin;
RP_result.U_RP_pxmax         = U_RP_pxmax;
RP_result.U_RP_pymin         = U_RP_pymin;
RP_result.U_RP_pymax         = U_RP_pymax;
RP_result.U_RP_pzbottom      = U_RP_pzbottom;



% ========================================================================
% Functions
% ========================================================================

function J_inv = calculate_inverse_jacobian( ...
    xixl, xiyl, xizl, ...
    etaxl, etayl, etazl, ...
    gammaxl, gammy, gammz, index)
% Build full 3D Jacobian matrix and invert it
% ------------------------------------------------------------------------
% ABOUT WHY THE INVERSE IS REQUIRED:
%
% In SPECFEM3D_Cartesian, the geometric derivatives stored in the binary
% mesh files (variables named xixl, xiyl, xizl, etc.) are *not* the
% standard Jacobian matrix entries ∂x/∂ξ, ∂x/∂η, ∂x/∂γ.
%
% Instead, SPECFEM stores the *inverse Jacobian*, meaning derivatives of
% the reference coordinates with respect to the physical coordinates:
%
%       xixl  = ∂ξ/∂x
%       xiyl  = ∂ξ/∂y
%       xizl  = ∂ξ/∂z
%       etaxl = ∂η/∂x
%       ...
%
% Therefore, the matrix built from SPECFEM fields is:
%
%       J_specfem = ∂(ξ, η, γ) / ∂(x, y, z)
%
% However, the representation-theorem and GLL surface integrals require
% the *standard* Jacobian:
%
%       J = ∂(x, y, z) / ∂(ξ, η, γ)
%
% which satisfies:
%
%       J = (J_specfem)^(-1)
%
% This is why we must explicitly compute the matrix inverse here — we are
% converting SPECFEM’s stored derivatives into the standard form required
% by the mathematical formulation of the RP boundary integral.
% ------------------------------------------------------------------------

    J = [ xixl(index),    xiyl(index),    xizl(index);
          etaxl(index),   etayl(index),   etazl(index);
          gammaxl(index), gammy(index),   gammz(index) ];

    J_inv = inv(J);
end


function J_inv = calculate_2d_jacobian_inverse(xixl, xiyl, etaxl, etayl, index)
% Build 2D Jacobian for XY-oriented surfaces

    J = [ xixl(index), xiyl(index);
          etaxl(index), etayl(index) ];

    J_inv = inv(J);
end


function J_inv = calculate_2d_jacobian_inverse_xz(xixl, xizl, gammaxl, gammazl, index)
% Build 2D Jacobian for XZ-oriented surfaces

    J = [ xixl(index), xizl(index);
          gammaxl(index), gammazl(index) ];

    J_inv = inv(J);
end


function J_inv = calculate_2d_jacobian_inverse_yz(etayl, etazl, gammay, gammaz, index)
% Build 2D Jacobian for YZ-oriented surfaces

    J = [ etayl(index), etazl(index);
          gammay(index), gammaz(index) ];

    J_inv = inv(J);
end


% ========================================================================
% Surface GLL Integral Computation
% ========================================================================
function [weights_matrix, area_integral] = compute_surface_integral( ...
    Ux_RP_all, ...
    xixl_all, xiyl_all, xizl_all, ...
    etaxl_all, etayl_all, etazl_all, ...
    gammaxl_all, gammayl_all, gammazl_all, ...
    wxgll_all, wygll_all, wzgll_all, ...
    point_indices, direction)
% Performs surface integration for a specific boundary face.
% direction = 'yz' | 'xz' | 'xy'

    [num_time_steps, num_space_points] = size(Ux_RP_all);

    area_integral = zeros(num_time_steps,1);
    weights_matrix = zeros(length(point_indices),1);

    for p = 1:length(point_indices)
        index = point_indices(p);

        % Choose 2D Jacobian form depending on face orientation
        if strcmp(direction,'yz')
            J_inv = calculate_2d_jacobian_inverse_yz( ...
                     etayl_all, etazl_all, ...
                     gammayl_all, gammazl_all, index );
            weight = wygll_all(index) * wzgll_all(index);

        elseif strcmp(direction,'xz')
            J_inv = calculate_2d_jacobian_inverse_xz( ...
                     xixl_all, xizl_all, ...
                     gammaxl_all, gammazl_all, index );
            weight = wxgll_all(index) * wzgll_all(index);

        elseif strcmp(direction,'xy')
            J_inv = calculate_2d_jacobian_inverse( ...
                     xixl_all, xiyl_all, ...
                     etaxl_all, etayl_all, index );
            weight = wxgll_all(index) * wygll_all(index);
        end

        % Handle invalid Jacobian values
        % NOTE:
        % This is not hard-coding. This correction is required to match the way
        % SPECFEM3D_Cartesian stores geometric derivatives for *regular elements*.
        %
        % In SPECFEM3D_Cartesian, for regular hexahedral elements the derivative
        % fields are stored in a compressed form:
        %
        %       ξ_x = (∂ξ/∂x)_regular,   but  ξ_y and ξ_z are NOT stored separately
        %
        % Because for regular elements:
        %
        %       ∂ξ/∂x = ξ_x_regular
        %       ∂ξ/∂y = 0
        %       ∂ξ/∂z = 0
        %
        % and similarly:
        %
        %       ∂η/∂x = 0,   ∂η/∂y = η_y_regular,   ∂η/∂z = 0
        %       ∂γ/∂x = 0,   ∂γ/∂y = 0,             ∂γ/∂z = γ_z_regular
        %
        % Because the SPECFEM output stores only the nonzero component, the other
        % two components must be explicitly set to zero. These assignments restore
        % the correct full derivative structure required for constructing the
        % geometric Jacobian.
        %
        if isnan(det(J_inv))
            J_det = 1/0.002^2; % must be changed when changing element sizes.
        else
            J_det = det(J_inv);
        end

        % Accumulate surface integral
        area_integral = area_integral + Ux_RP_all(:,p) * J_det * weight;
        weights_matrix(p) = J_det * weight;
    end

end


end % End of main function
