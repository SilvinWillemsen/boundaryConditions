close all;
testPrev = true;
if ~testPrev
    clear all;
    testPrev = false;
end
Nold = 14;
Nnew = 10;

if testPrev
    NoldTmp = Nold;
    Nold = Nnew;
    Nnew = NoldTmp;
end
alpha = 0;
bound = "none";
interpolType = "cubic";
% Nold = 10;
% Nnew = 19;

% If spatial upsampling, don't do anything fancy
if Nold < Nnew
    flag = false;
else % otherwise
    flag = true;
    NoldTmp = Nold;
    Nold = Nnew;
    Nnew = NoldTmp;
end

alphaInc = Nold / Nnew;
if bound == "clamped"
    alpha = alpha + 2 * alphaInc;
    iRange = 3:Nnew-2;
elseif bound == "ss"
    alpha = alpha + alphaInc;
    iRange = 2:Nnew-1;
else
    iRange = 1:Nnew;
end

I = zeros(Nnew, Nold);
j = 1;
for i = iRange
    if alpha >= 1
        j = j + floor(alpha);
        alpha = mod(alpha, 1);
    end
    if interpolType == "linear"
        interpolVect = [(1-alpha); alpha];
        if j >= Nold
            I(i, Nold) = interpolVect(1);
        else
            I(i, j:j+1) = interpolVect;
        end
    elseif interpolType == "cubic"
        interpolVect = [(alpha * (alpha - 1) * (alpha - 2)) / -6.0; ...
            ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2.0; ...
            (alpha * (alpha + 1) * (alpha - 2)) / -2.0; ...
            (alpha * (alpha + 1) * (alpha - 1)) / 6.0];
        if j == 1
        I(i, j:j+2) = interpolVect(2:end);
        elseif j >= Nold-1
            I(i,j-1:Nold) = interpolVect(1:(Nold-(j-2)));
        else
            I(i, j-1:j+2) = interpolVect;
        end
    end

    alpha = alpha + alphaInc;
end
if flag
    Nnew = Nold;
    Nold = NoldTmp;
    Iinv = I;
    I = I' * Nnew/Nold;
end

u = sin(2 * pi * (0:Nold-1)' / Nold);

plot((0:Nold-1)/Nold, u);
uInter = I * u;
hold on;
if testPrev
    plot((0:Nnew-1)/Nnew, uInter);
    plot((0:Nold-1)/Nold, uInterSave);
    legend(["Original", "Interpolated", "Reinterpolated"])
else
    uInterSave = uInter;
end

% end