addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");
addpath("../exercises");

%% t1
close all

N = 10000;
data = randn(1, N);

A = [1 -1.79 0.84]; % AR (on y)
C = [1 -0.81 -0.11]; % MA (on e)

figure
plotACFnPACFnReal(A, 1, data, 1, "AR")
plotACFnPACFnReal(1, C, data, 2, "MA")
plotACFnPACFnReal(A, C, data, 3, "ARMA")

% With the plots we can identify an AR, MA and ARMA/AR processes.
% When N is changed to 100 it harder to identify the process.

%% t2

close all

N = 10000;
data = randn(1, N);

A = [1 -1.79 0.84]; % AR (on y)
C = [1 -0.81 -0.11]; % MA (on e)

Aoutside = [1 -2 0.84];
Coutside = [1 -1 -0.11];

figure
plotRealnOustide(A, C, Aoutside, Coutside, data);

% cannot really see any difference for C outside.

%% t3

close all

load week2data
data = y;

% figure
% plot(data)
% ACFnPACFnNormplot(data, 32);

A = [1 0 0 0]; B = []; C = [1 0 0];
model_init = idpoly(A, B, C);
model_arma = pem(data, model_init);
figure
resid(model_arma, data);
rarma = resid(model_arma, data);
ACFnPACFnNormplot(rarma.OutputData, 32);

checkIfNormal(rarma.OutputData, "model_arma");
checkIfWhite(rarma.OutputData);

%% functions

function plotACFnPACFnReal(A, C, data, plotInd, procStr)
    noLags = 32;
    signLvl = 0.05;

    dataFiltered = filter(C, A, data);
    dataFiltered = dataFiltered(10:end);

    subplot(3, 3, plotInd)
    plot(dataFiltered);
    title(sprintf("Data realization of a %s-process", procStr))

    subplot(3, 3, plotInd + 3)
    acfEst = acf(dataFiltered, noLags, signLvl, 1);
    title(sprintf("ACF %s", procStr))

    subplot(3, 3, plotInd + 6)
    pacfEst = pacf(dataFiltered, noLags, signLvl, 1);
    title(sprintf("PACF %s", procStr))
end

function plotRealnOustide(A, C, Aoutside, Coutside, data)
    ar = filter(1, A, data);
    arOutside = filter(1, Aoutside, data);

    subplot(4, 3, 1); plot(ar); title("AR");
    subplot(4, 3, 4); plot(arOutside); title("AR outside");

    ma = filter(C, 1, data);
    maOutside = filter(Coutside, 1, data);

    subplot(4, 3, 2); plot(ma); title("MA");
    subplot(4, 3, 5); plot(maOutside); title("MA outside");

    arma = filter(C, A, data);
    armaAoutside = filter(C, Aoutside, data);
    armaCoutside = filter(Coutside, A, data);
    armaOutside = filter(Coutside, Aoutside, data);

    subplot(4, 3, 3); plot(arma); title("ARMA");
    subplot(4, 3, 6); plot(armaAoutside); title("ARMA A outside");
    subplot(4, 3, 9); plot(armaCoutside); title("ARMA C outside");
    subplot(4, 3, 12); plot(armaOutside); title("ARMA outside");
end
