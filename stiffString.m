clear all;
close all;

exc = "cos";

fs = 44100;
k = 1/fs;

f0 = 196;
lengthSound = fs;

%% String parameters
L = 1;
r = 0.0005;
rho = 7850;
A = r^2 * pi;

% Wavespeed
c = f0 * 2;
T = c * c * rho * A; 

% Stiffness
E = 2e11;
I = r^4 * pi / 4;

% Damping
s0 = 0.1;
s1 = 0.005;

% Boundary
bound = "ss";

kappa = sqrt(E * I / (rho * A));

[B, C, N, h, Dxx, Dxxxx, s0, s1] = unscaledStringBoundaryCond(rho, A, T, E, I, L, s0, s1, k, bound);

%% State vectors

uNext = zeros(N,1);
u = zeros(N,1);

if exc == "cos"
    loc = floor(0.5 * N);
    width = floor(floor(N / 10) / 2) * 2;
    raisedCos = (1 - cos(2 * pi * [0:width] / width)) * 0.001;
    startIdx = loc - width / 2;
    endIdx = loc + width / 2;
    u(startIdx:endIdx) = raisedCos;
end
uPrev = u;

%% Energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
energy = zeros(lengthSound, 1);

eVec = 2:N-1;
eVecTest = 1:N-1;
vec = 2:N-1;

interpolate = false;

for n = 1:lengthSound
    uNext = B * u + C * uPrev;
    
    %% Energystuff
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * 1/h * sum((u(eVecTest+1) - u(eVecTest)) .* (uPrev(eVecTest+1) - uPrev(eVecTest)))...
        + E * I / 2 * 1/h^3 * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
    energy(n) = kinEnergy(n) + potEnergy(n);
    
    rOCkinEnergy(n) = h * rho * A / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = h * T / (2*k*h^2) * sum((Dxx * u) .* (uNext - uPrev)) ...
         - h * E * I / (2 * k * h^4) * sum(Dxxxx * u .* (uNext - uPrev));%...
    rOCdamp0StringEnergy(n) = -2 * s0 * h / (4 * k^2) * sum((uNext - uPrev).*(uNext - uPrev));
    rOCdamp1StringEnergy(n) = 2 * h * s1 / (2 * k^2 * h^2) * sum((Dxx * u - Dxx * uPrev) .* (uNext - uPrev));
    rOCenergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdamp0StringEnergy(n) - rOCdamp1StringEnergy(n);
    
    if n = 5000
        interpolate = true;
    end
    
    if interpolate
        [B, C, Ntmp, h, Dxx, Dxxxx, s0, s1] = unscaledStringBoundaryCond(rho, A, T, E, I, L, s0, s1, k, bound);
        I = 
        interpolate = false;
    end
    %% Update states
    uPrev = u;
    u = uNext;
    
    subplot(2,1,1);
    plot(uNext);
    subplot(2,1,2)
%     plot(energy(1:n)/energy(1) - 1);
    plot(rOCenergy(1:n));
    drawnow;
end


