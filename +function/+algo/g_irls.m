function [estimated_pt, iter] = g_irls(b,A,lambda, maxiter, p, tolerance, s_amount,vs_xyz, min_d)
% Function  : g-irls
% Creation  : 01-11-25
% Author    : Weiting Lai
% Version   : v1.0 (01-11-25)
%
% Description
%
%   Implement Group Iteratively Reweighted Least Squares (G-IRLS) from the 
%   reference paper.
%
% Reference 
% 
% S. Koyama et al., Sparse sound field decomposition for super-resolution 
% in recording and reproduction, J. Acoust. Soc. Amer., 2018.

  % Initialize the solution x
  N = size(A,2);
  F = size(A,3);
  x = zeros(N,1,F);
  eta = vecnorm(x, 2, 3);
  w = (eta.^2+lambda).^(p/2-1);
  
  Q = spdiags(1./w(:), 0, N, N);
  
  % Iterate until convergence
    for i = 1 : maxiter
        
        x_hat = x;
        x_temp = x_hat;
        % Solve the weighted least squares problem
        for f = 1:F
            temp = pinv(A(:,:,f)*Q*A(:,:,f)');
            x_hat(:,:,f) = Q*A(:,:,f)'*temp*b(:,:,f);
        end
        eta = vecnorm(x_hat, 2, 3);
        % Compute the weight vector w
        w = (eta.^2+lambda).^(p/2-1);
        Q = diag(1./w);

        % Check for convergence
        tol = norm(x_hat(:) - x_temp(:), 1);
        if tol < tolerance
            break
        end
    end
    iter = i;
    x_hat = abs(sum(x_hat,3));
    % sort the x_hat
    [~, index] = sort(x_hat,'descend');

    % initialisation min_d & est_pt
    pti = zeros(s_amount,1);
    estimated_pt = zeros(s_amount,3);
    pti(1) = index(1);
    estimated_pt(1,:) = vs_xyz(pti(1),:);

    % remove points shorter than min_d
    left = 1;
    right = 2;
    while left < s_amount
        if_allowed = 1;
        for i = 1:left
            % distance between found points and the next largest point
            dist = norm(estimated_pt(i,:)-vs_xyz(index(right),:));
            if(dist < min_d)
                if_allowed = 0;
                break;
            end
        end
        % store new point
        if(if_allowed)
            left = left + 1;
            pti(left) = index(right);
            estimated_pt(left,:) = vs_xyz(pti(left),:);
        end
        right = right + 1;
    end
end