function [B, C, N, h, Dxx, Dxxxx, s0, s1] = unscaledStringBoundaryCond(rho, A, T, E, I, L, s0, s1, k, bound)
   
    kappa = sqrt(E * I / (rho * A));
    c = sqrt(T/(rho * A));

    h = sqrt((c^2 * k^2 + 4 * s1 * k + sqrt((c^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
    N = floor(L/h); % Number of gridpoints
    h = L/N; % Recalculate gridspacing
    
    s0 = s0 * rho * A;
    s1 = s1 * rho * A;
    Ndiff = N;
    if bound == "clamped"
        N = N - 4;
        midVect = 6 * ones(1, N);
    elseif bound == "ss"
        N = N - 2;
        midVect = [5, 6 * ones(1, N-2), 5];
    end
    
    Ndiff = Ndiff - N;

    Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, midVect, N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
    Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));

    B = (2 * rho * A / k^2 * eye(N) + T / h^2 * Dxx - E * I / h^4 * Dxxxx + 2 * s1 / (k * h^2) * Dxx) / (rho * A / k^2 + s0 / k);
    C = ((-rho * A / k^2 + s0 / k) * eye(N) - 2 * s1 / (k*h^2) * Dxx) / (rho * A / k^2 + s0/k);

    N = N + Ndiff;
    
    halfNdiff = Ndiff / 2;
    Dxx2 = zeros(N);
    Dxx2(1+halfNdiff:end-halfNdiff, 1+halfNdiff:end-halfNdiff) = Dxx;
    Dxx = Dxx2;
    
    Dxxxx2 = zeros(N);
    Dxxxx2(1+halfNdiff:end-halfNdiff, 1+halfNdiff:end-halfNdiff) = Dxxxx;
    Dxxxx = Dxxxx2;
    
    B2 = zeros(N);
    B2(1+halfNdiff:end-halfNdiff, 1+halfNdiff:end-halfNdiff) = B;
    B = B2;

    C2 = zeros(N);
    C2(1+halfNdiff:end-halfNdiff, 1+halfNdiff:end-halfNdiff) = C;
    C = C2;
end