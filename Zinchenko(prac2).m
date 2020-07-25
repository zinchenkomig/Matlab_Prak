%% 0
T = 8;
k = 3;
M = 3;
L = 5;
eps = 10;
NGrid = 1000;
BeginN = 2;
switchGrid = 50;
deltaEndSet = 0.3;

x1MinMatr = zeros(NGrid, 10);
x2MinMatr = zeros(NGrid, 10);
psi1MinMatr = zeros(NGrid, 10);
psi2MinMatr = zeros(NGrid, 10);
u1MinMatr = zeros(NGrid, 10);
u2MinMatr = zeros(NGrid, 10);
SwitchTime = zeros(2, 10);
%% 1
% psi1 > 0; psi2 >= 0
% Two switches
tVec = linspace(0, T, switchGrid);
tswitch1 = tVec(2:end);
tVec = tswitch1;

A = repelem(tVec, 1, size(tVec, 2))';
B = repmat(tVec, 1, size(tVec, 2))';
C = A < B;
tswitch1 = A(C);
tswitch2 = B(C);

swSize = size(tswitch1,1);

psi1 = 4*k^3./((exp(k*(tswitch2 - tswitch1)) - 1).*(exp(k*tswitch1) - 1).^2);
%all possible psi1 for switch moments


t = linspace(0,T,NGrid)';
t = repmat(t,1,size(tswitch1,1))';

psi21 = (psi1/k).*(exp(k*(tswitch1 - t)) - 1);



psi22 = (psi1/k).*((exp(k*(tswitch1 - tswitch2)) - 2).*exp(k*(t - tswitch2)) + 1);
psi21(t>tswitch2) = psi22 (t> tswitch2);


x21 = (psi1/(4*k^2)) .* (2 - 2*exp(k*t) - exp(k*(tswitch1 - t)) + exp(k*(t + tswitch1)));
x22 = exp(k*(t-tswitch1)) .* (psi1/(4*k^2)) .* (exp(k*tswitch1) - 1).^2;
x23 = (psi1/(4*k^2)) .* exp(k*(2*tswitch2 - tswitch1 - t)) .* (exp(k*tswitch1) - 1).^2;
x21(t>tswitch1 & t<tswitch2) = x22(t>tswitch1 & t<tswitch2);
x21(t>tswitch2) = x23(t>tswitch2);

x10 = linspace(-M, M, BeginN);
x11 = zeros(swSize, NGrid, BeginN);

x1t1 = zeros(swSize,2);
x1t2 = zeros(swSize,2);

for i = 1:BeginN
    xbuf1 = x10(i) + (psi1/(4*k^3)).*(2*k*t - 2*exp(k*t) + 2 + exp(k*(tswitch1 - t)) - 2*exp(k*tswitch1) + exp(k*(t+tswitch1)))+ k*t;
    x1t1(:,i) = x10(i) + (psi1/(4*k^3)).*(2*k*tswitch1 - 4*exp(k*tswitch1) + 3 + exp(2*k*(tswitch1)))+ k*tswitch1;
    
    xbuf2 = x1t1(:,i) + ...
        (psi1/(4*k^3)) .* (exp(k*tswitch1) - 1).^2 .*( exp(k*(t - tswitch1)) - 1) + ...
        k*(t-tswitch1);
    x1t2(:,i) = x1t1(:,i) + (psi1/(4*k^3)) .* (exp(k*tswitch1) - 1).^2 .*( exp(k*(tswitch2 - tswitch1)) - 1) + ...
        k*(tswitch2-tswitch1);
    xbuf3 = x1t2(:,i) + ...
        (psi1/(4*k^3)) .* (exp(k*tswitch1) - 1).^2 .* (exp(k*(tswitch2-tswitch1)) - ...
        exp(k*(2*tswitch2 - tswitch1 - t))) - k*(t-tswitch2);
    xbuf1(t>tswitch1 & t<tswitch2) = xbuf2(t>tswitch1 & t<tswitch2);
    xbuf1(t>tswitch2) = xbuf3(t>tswitch2);
    x11(:,:,i) = xbuf1;
end



x1Matr = x11(:, :, 1)';
x2Matr = x21(:, :)';

psi2Matr = psi21';
psi1Matr = repmat(psi1', NGrid, 1);
psi21(t>tswitch1) = 0;
psi21 = psi21';

u1Matr = psi21 / 2;
u2Matr = zeros(size(u1Matr));
u2Matr(t'<tswitch2') = k;
u2Matr(t'>=tswitch2') = -k;

x1Matr = [x1Matr, x11(:, :, 2)'];
x2Matr = [x2Matr, x21(:, :)'];
IndVec = find(abs(x1Matr(end,:)-L) < deltaEndSet & abs(x2Matr(end,:))<eps);


psi1Matr = [psi1Matr, psi1Matr];
psi2Matr = [psi2Matr, psi2Matr];
tswitch1 = [tswitch1; tswitch1];
tswitch2 = [tswitch2; tswitch2];
u1Matr = [u1Matr, u1Matr];
u2Matr = [u2Matr, u2Matr];

J1 = NaN;
if ~isempty(IndVec)
[J1, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);

r = 1;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:,minN);
u2MinMatr(:,r) = u2Matr(:, minN);
SwitchTime(:,r) = [tswitch1(minN), tswitch2(minN)];
end
end
%% 2
% psi1 > 0; psi2 <= 0

% No switches

t = linspace(0,T,NGrid);
x1Matr = [-M + k*t', M + k*t'];
x2Matr = zeros(size(x1));
IndVec = find(abs(x1Matr(end,:)-L) < deltaEndSet);
r = 2;
J2 = NaN;
if ~isempty(IndVec)
    J2 = 0;
    x1MinMatr(:,r) = x1Matr(:,IndVec(1));
    x2MinMatr(:,r) = x2Matr(:,IndVec(1));
    psi1MinMatr(:,r) = 0;
    psi2MinMatr(:,r) = 0;
    u1MinMatr(:,r) = 0;
    u2MinMatr(:,r) = k;
    SwitchTime(:,r) = [NaN, NaN];
end


%% 3
tVec = linspace(0, T, switchGrid);
tswitch1 = tVec(2:end);
tVec = tswitch1;

A = repelem(tVec, 1, size(tVec, 2))';
B = repmat(tVec, 1, size(tVec, 2))';
C = A < B;
tswitch1 = A(C);
tswitch2 = B(C);

swSize = size(tswitch1,1);

psi1 = 4*k^3 ./ ((2 - exp(k*(tswitch1 - tswitch2)) - exp(k*(tswitch2 - tswitch1))) .* (1 - exp(k*(tswitch2 - tswitch1))));
%all possible psi1 for switch moments


t = linspace(0,T,NGrid)';
t = repmat(t,1,size(tswitch1,1))';

psi21 = (psi1/k).*(1 - exp(k*(t - tswitch1)));


psi2t2 = (psi1/k) .* (1 - exp(k*(tswitch2 - tswitch1)));
psi22 = exp(k*(tswitch2 - t)) .* (psi2t2 + (psi1/k)) - psi1/k;
psi21(t>=tswitch2) = psi22 (t>=tswitch2);


x21 = zeros(size(t));
x22 = (psi1/(4*k^2)) .* (2 - exp(k*(tswitch1-t)) - exp(k*(t - tswitch1)));
x2t2 = (psi1/(4*k^2)) .* (2 - exp(k*(tswitch1-tswitch2)) - exp(k*(tswitch2 - tswitch1)));

x23 = x2t2 .* exp(k*(t - tswitch2)) + (1/(2*k)) * sinh(k*(t-tswitch2)) .* (psi2t2 + psi1/k) +...
    (psi1/(2*k^2)) .* (1 - exp(k*(t - tswitch2)));

x21(t>tswitch1 & t<tswitch2) = x22(t>tswitch1 & t<tswitch2);
x21(t>=tswitch2) = x23(t>=tswitch2);

x10 = linspace(-M, M, BeginN);
x11 = zeros(swSize, NGrid, BeginN);

x1t1 = zeros(swSize,2);
x1t2 = zeros(swSize,2);
for i = 1:BeginN
    xbuf1 = x10(i) - k*t;
    x1t1(:,i) = x10(i)-k*tswitch1;

    xbuf2 = x1t1(:,i) + ...
        (psi1/(4*k^2)) .* (2*(t-tswitch1) + (2/k) * sinh(k*(tswitch1- t))) - k*(t-tswitch1);

    x1t2(:,i) = x1t1(:,i) + ...
        (psi1/(4*k^2)) .* (2*(tswitch2-tswitch1) + (2/k) * sinh(k*(tswitch1- tswitch2))) - k*(tswitch2-tswitch1);

    xbuf3 = x1t2(:,i) + ...
       (x2t2/k).* (exp(k*(t - tswitch2)) - 1) + ...
       (1/(2*k^2)) * (psi2t2 + psi1/k).*(cosh(k*(t-tswitch2))  - 1) +...
       (psi1/(2 * k^2)) .* (t-tswitch2) -...
       (psi1/(2*k^3)) .* (exp(k*(t - tswitch1)) - exp(k*(tswitch2 - tswitch1)))+...
       k * (t - tswitch2);

   xbuf1(t>tswitch1 & t<tswitch2) = xbuf2(t>tswitch1 & t<tswitch2);
   xbuf1(t>=tswitch2) = xbuf3(t>=tswitch2);
   x11(:,:,i) = xbuf1;
end

x1Matr = x11(:, :, 1)';
x2Matr = x21(:, :)';

psi2Matr = psi21';
psi1Matr = repmat(psi1', NGrid, 1);
psi21(t<tswitch1) = 0;
psi21 = psi21';

u1Matr = psi21 / 2;
u2Matr = zeros(size(u1Matr));
u2Matr(t'<tswitch2') = -k;
u2Matr(t'>=tswitch2') = k;

x1Matr = [x1Matr, x11(:, :, 2)'];
x2Matr = [x2Matr, x21(:, :)'];

IndVec = find(abs(x1Matr(end,:)-L) < deltaEndSet & abs(x2Matr(end,:))<eps);

psi1Matr = [psi1Matr, psi1Matr];
psi2Matr = [psi2Matr, psi2Matr];
tswitch1 = [tswitch1; tswitch1];
tswitch2 = [tswitch2; tswitch2];
u1Matr = [u1Matr, u1Matr];
u2Matr = [u2Matr, u2Matr];

J3 = NaN;
if ~isempty(IndVec)
[J3, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);

r = 3;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:,minN);
u2MinMatr(:,r) = u2Matr(:, minN);
SwitchTime(:,r) = [tswitch1(minN), tswitch2(minN)];
end
end

%% 4

% psi1 < 0; psi2 > 0

tVec = linspace(0, T, switchGrid);
t = linspace(0, T, NGrid);

x1Matr = zeros(NGrid, switchGrid);
x2Matr = zeros(NGrid, switchGrid);
psi2Matr = zeros(NGrid, switchGrid);
u2Matr = zeros(size(psi2Matr));
x10 = M;
x2t2 = eps;
psi1Vec = zeros(1, switchGrid);

for j = 2:switchGrid
    
    tswitch1 = tVec(j);
    
    
    Fpsi2t1 = @(psi1) (psi1/k)*(1 - exp(-2*k*tswitch1) - 2*exp(-k*tswitch1) + (2*k^2/psi1) * x2t2 * exp(k*(tswitch1 - T)) + exp(2*k*(tswitch1 - T)) - exp(2*k*tswitch1) - 2*exp(k*(tswitch1 - T)))/(1 - exp(-2*k*tswitch1) - exp(2*k*(tswitch1 - T)) + exp(2*k*tswitch1));
    Fx2t1 = @(psi1) (1/(2*k)) * (Fpsi2t1(psi1) - psi1/k)*(1 - exp(-2*k*tswitch1)) + (psi1/(k^2)) * (1 - exp(-k*tswitch1));
    psi1 = fsolve(@(psi1) psi1 + Fpsi2t1(psi1) * Fx2t1(psi1), 1);
    
    psi2t1 = Fpsi2t1(psi1);
    x21 = (psi2t1 - psi1/k)*(1/(2*k)) * (exp(k*(t-tswitch1)) - exp(-k*(t+tswitch1))) + (psi1/(k^2))*(1 - exp(-k*t));
    x2t1 = Fx2t1(psi1);
    x22 = x2t2 * exp(k*(t - T)) + (1/(2*k)) * (psi2t1 + psi1/k) * (exp(k*(t + tswitch1 - 2*T)) - exp(k*(t + tswitch1))) + (psi1/(k^2)) * (1 - exp(k*(t - T)));
    x21(t>=tswitch1) = x22(t>=tswitch1);
    
    psi21 = exp(k*(t - tswitch1)) * psi2t1 + (psi1/k)*(1 - exp(k*(t - tswitch1)));
    psi22 = exp(-k*(t - tswitch1)) * psi2t1 + (psi1/k) * (exp(k*(tswitch1 - t)) - 1);
    psi21(t>=tswitch1) = psi22(t>=tswitch1);
    
    tVec1 = t(t<tswitch1);
    
    beforeSw = size(tVec1,2);
    x11 = zeros(size(tVec1));
    for i = 2:beforeSw
        x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
    end
    x11(1) = x11(2);
    
    tVec2 = t(t>=tswitch1);
    
    x12 = zeros(size(tVec2));
    for i = 2:size(tVec2,2)
        x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
    end
    x12(1) = x11(end);
    
    x1 = [x11, x12];
    
    x1Matr(:,j) = x1;
    x2Matr(:,j) = x21;
    psi2Matr(:,j) = psi21;
    psi1Vec(j) = psi1;
    u2Matr(t<tswitch1,j) = -k;
    u2Matr(t>=tswitch1,j) = k;
end

x1Matr = [x1Matr, x1Matr - 2*M];
x2Matr = [x2Matr, x2Matr];
psi2Matr = [psi2Matr, psi2Matr];
u1Matr = psi2Matr/2;

psi1Matr = repmat(psi1Vec, NGrid, 1);

IndVec = find(abs(x1Matr(end,:,1)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);

J4 = NaN;
if ~isempty(IndVec)
    
[J4, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 4;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:,minN);
u2MinMatr(:,r) = u2Matr(:, minN);
SwitchTime(:,r) = [tVec(minN), NaN];
end

end



%% 5

% psi1 = 0; psi2 > 0
t = linspace(0, T, NGrid);

psi20 = 2*eps*k/(sinh(k*T));
x10 = M;

x2 = (psi20/(2*k)) * sinh(k*t);
psi2 = psi20 * exp(-k*t);
psi1 = zeros(size(t));
x1 = x10*((psi20/(2*k^2)) * (cosh(k*t) - 1) + k*t);

r = 5;

if (abs(x1(end)-L) < deltaEndSet && abs(x2(end))<=eps)
    u1Matr = psi2/2;
    J5 = trapz(t, u1Matr.^2);
    x1MinMatr(:, r) = x1;
    x2MinMatr(:, r) = x2;
    psi1MinMatr(:, r) = psi1;
    psi2MinMatr(:, r) = psi2;
    u1MinMatr(:,r) = 0;
    u2MinMatr(:,r) = k;
    SwitchTime(:,r) = [NaN, NaN];
else
    J5 = NaN;
end


%% second u


tSwitches = linspace(0,T,switchGrid);
t = linspace(0,T,NGrid);

x1Matr = zeros(NGrid, 4*switchGrid);
x2Matr = zeros(NGrid, 4*switchGrid);
psi1Matr = zeros(NGrid, 4*switchGrid);
psi2Matr = zeros(NGrid, 4*switchGrid);
u2Matr = zeros(NGrid, 4*switchGrid);
% Transvers: x2(T) == +-eps; psi2(T) != 0;
edgesMatr =...
    [M, eps;...
    M, -eps;...
    -M, eps;...
    -M, -eps];
for l = 1:4
    
    x10 = edgesMatr(l,1);
    x2T = edgesMatr(l,2);
    for j = 2:switchGrid
        
        tswitch1 = tSwitches(j);
        
        psi = fsolve(@(psi) psi2x2Sys1(psi, tswitch1,x2T,T,k), [1,1]);
        if psi(1) <= 0
            continue;
        end
        x21 =(1/(4*k))*(psi(2) + (psi(1)/k))*(exp(k*(t+tswitch1)) - exp(k*(tswitch1 - t))) + (psi(1)/(2*k^2)) * (1 - exp(k * t)); 
        x22 = x2T*exp(-k*(t - T)) + (1/(4*k)) * (psi(2) - (psi(1)/k)) * (exp(k*(t - tswitch1)) - exp(k*(2*T - t - tswitch1))) - (psi(1)/(2*k^2)) * (exp(-k*(t - T)) - 1);
        
        psi21 = psi(2) * exp(-k*(t - tswitch1)) + (psi(1)/k) * (exp(k*(tswitch1 - t)) - 1);
        psi22 = psi(2) * exp(k*(t - tswitch1)) + (psi(1)/k) * (1 - exp(k*(t - tswitch1)));
        
        x21(t>=tswitch1) = x22(t>=tswitch1);
        psi21(t>=tswitch1) = psi22(t>=tswitch1);
        
        
        tVec1 = t(t<tswitch1);
        
        beforeSw = size(tVec1,2);
        x11 = zeros(size(tVec1));
        for i = 2:beforeSw
            x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
        end
        x11(1) = x11(2);
        
        tVec2 = t(t>=tswitch1);
        
        x12 = zeros(size(tVec2));
        for i = 2:size(tVec2,2)
            x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
        end
        x12(1) = x11(end);
        if size(x12, 2) > 1
            x12(1) = (x11(end) + x12(2)) / 2;
        end
        
        x1 = [x11, x12];
        u2 = zeros(1, NGrid);
        u2(t<tswitch1) = k;
        u2(t>=tswitch1) = -k;
        x1Matr(:,j + (l-1)*switchGrid) = x1;
        x2Matr(:,j+ (l-1)*switchGrid) = x21;
        psi2Matr(:,j+ (l - 1)*switchGrid) = psi21;
        psi1Matr(:, j+ (l - 1)*switchGrid) = psi(1);
        u2Matr(:,j + (l-1)*switchGrid) = u2;
    end
end

u1Matr = psi2Matr / 2;
IndVec = find(abs(x1Matr(end,:,1)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);
J6 = NaN;

if ~isempty(IndVec)
    
[J6, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 6;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:, minN);
u2MinMatr(:,r) = u2Matr(:, minN);
tSwitches = repmat(tSwitches, 1, 4);
SwitchTime(:,r) = [tSwitches(minN), NaN];
end

end



%% Special 1

t = linspace(0,T,NGrid);
syms t1;
x2T = -eps;
        psi(2) = -2*eps^2;
        psi(1) = x2T * psi(2);
        tswitch1 = solve((1/(4*k))*((2*eps^3)/k - 2*eps^2)*(exp(2*k*t1) - 1) + ((eps^3)/(k^2)) * (1 - exp(k * t1)) + eps == 0, t1);
        tswitch1 = double(tswitch1);
        ts1 = tswitch1(tswitch1>=0 & tswitch1<T);
        
        if ~isempty(ts1)
        ts1 = ts1(isreal(ts1));
        end
        
switchAm = size(ts1,2);
x1Matr = zeros(NGrid, 2*switchAm);
x2Matr = zeros(NGrid, 2*switchAm);
psi1Matr = zeros(NGrid, 2*switchAm);
psi2Matr = zeros(NGrid, 2*switchAm);

J7 = NaN;

edgesMatr =...
    [-M;M];
x2T = -eps;
if ~isempty(ts1)
for l = 1:2
    
    x10 = edgesMatr(l);
    
        for j = 1:size(t1,2)
            tswitch1 = ts1(j);
        x21 = (1/(4*k))*((2*eps^3)/k - 2*eps^2)*(exp(k*(t + tswitch1)) - exp(k*(tswitch1 - t))) + ((2*eps^3)/(2*k^2)) * (1 - exp(k * t));
        psi21 = psi(2) * exp(-k*(t - tswitch1)) + (psi(1)/k) * (exp(k*(tswitch1 - t)) - 1);
        
        x21(t>=tswitch1) = x2T;
        psi21(t>=tswitch1) = psi(2);
        
        
        tVec1 = t(t<tswitch1);
        
        beforeSw = size(tVec1,2);
        x11 = zeros(size(tVec1));
        for i = 2:beforeSw
            x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
        end
        if (size(x11,2) ~= 1)
            x11(1) = x11(2);
        else
            x11(1) = x12(1);
        end
        tVec2 = t(t>=tswitch1);
        
        x12 = zeros(size(tVec2));
        for i = 2:size(tVec2,2)
            x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
        end
        x12(1) = x11(end);
        
        x1 = [x11, x12];
        u2 = zeros(1, NGrid);
        u2(t<tswitch1) = k;
        u2(t>=tswitch1) = x2T;
        
        x1Matr(:,j + (l-1)) = x1;
        x2Matr(:,j+ (l-1)) = x21;
        psi2Matr(:,j+ (l - 1)) = psi21;
        psi1Matr(:, j+ (l - 1)) = psi(1);
        u2Matr(:,j + (l-1)*switchGrid) = u2;
        end
end


u1Matr = psi2Matr / 2;
IndVec = find(abs(x1Matr(end,:)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);


if ~isempty(IndVec)
    
[J7, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 7;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:, minN);
u2MinMatr(:,r) = u2Matr(:, minN);
tSwitches = repmat(tSwitches, 1, 4);
SwitchTime(:,r) = [tSwitches(minN), NaN];
end

end

end

%% psi1 < 0

tSwitches = linspace(0,T,switchGrid);
t = linspace(0,T,NGrid);

x1Matr = zeros(NGrid, 4*switchGrid);
x2Matr = zeros(NGrid, 4*switchGrid);
psi1Matr = zeros(NGrid, 4*switchGrid);
psi2Matr = zeros(NGrid, 4*switchGrid);
u2Matr = zeros(NGrid, 4*switchGrid);

% Transvers: x2(T) == +-eps; psi2(T) != 0;
edgesMatr =...
    [M, eps;...
    M, -eps;...
    -M, eps;...
    -M, -eps];
for l = 1:4
    
    x10 = edgesMatr(l,1);
    x2T = edgesMatr(l,2);
    for j = 2:switchGrid
        
        tswitch1 = tSwitches(j);
        
        psi = fsolve(@(psi) psi2x2Sys2(psi, tswitch1,x2T,T,k), [1,1]);
        if psi(1) > 0
            continue;
        end
        x21 =(1/(4*k))*(psi(2) - (psi(1)/k))*(exp(k*(t-tswitch1)) - exp(-k*(tswitch1 + t))) - (psi(1)/(2*k^2)) * (exp(-k * t) - 1); 
        x22 = x2T*exp(k*(t - T)) + (1/(4*k)) * (psi(2) + (psi(1)/k)) * (exp(k*(t + tswitch1 - 2*T)) - exp(k*(tswitch1 - t))) + (psi(1)/(2*k^2)) * (1 - exp(k*(t - T)));
        
        psi21 = psi(2) * exp(k*(t - tswitch1)) + (psi(1)/k) * (1 - exp(k*(t - tswitch1)));
        psi22= psi(2) * exp(-k*(t - tswitch1)) + (psi(1)/k) * (exp(k*(tswitch1 - t)) - 1);
        
        x21(t>=tswitch1) = x22(t>=tswitch1);
        psi21(t>=tswitch1) = psi22(t>=tswitch1);
        
        
        tVec1 = t(t<tswitch1);
        
        beforeSw = size(tVec1,2);
        x11 = zeros(size(tVec1));
        for i = 2:beforeSw
            x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
        end
        x11(1) = x11(2);
        
        tVec2 = t(t>=tswitch1);
        
        x12 = zeros(size(tVec2));
        for i = 2:size(tVec2,2)
            x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
        end
        x12(1) = x11(end);
        if size(x12, 2) > 1
            x12(1) = (x11(end) + x12(2)) / 2;
        end
        
        x1 = [x11, x12];
        
        u2 = zeros(1, NGrid);
        u2(t<tswitch1) = -k;
        u2(t>=tswitch1) = k;
        x1Matr(:,j + (l-1)*switchGrid) = x1;
        x2Matr(:,j+ (l-1)*switchGrid) = x21;
        psi2Matr(:,j+ (l - 1)*switchGrid) = psi21;
        psi1Matr(:, j+ (l - 1)*switchGrid) = psi(1);
        u2Matr(:,j + (l-1)*switchGrid) = u2;
    end
end


u1Matr = psi2Matr / 2;
IndVec = find(abs(x1Matr(end,:,1)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);

J8 = NaN;

if ~isempty(IndVec)
    
[J8, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 8;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:, minN);
u2MinMatr(:,r) = u2Matr(:, minN);
tSwitches = repmat(tSwitches, 1, 4);
SwitchTime(:,r) = [tSwitches(minN), NaN];
end

end


%% trans2


tSwitches = linspace(0,T,switchGrid);
t = linspace(0,T,NGrid);

x1Matr = zeros(NGrid, 2*switchGrid);
x2Matr = zeros(NGrid, 2*switchGrid);
psi1Matr = zeros(NGrid, 2*switchGrid);
psi2Matr = zeros(NGrid, 2*switchGrid);
u2Matr = zeros(NGrid, 2*switchGrid);

 %Transvers: psi2(T) = 0; x2(T) != +-eps;
edgesMatr = [-M; M];

for l = 1:2
    
    x10 = edgesMatr(l);
    
    for j = 2:switchGrid
        
        tswitch1 = tSwitches(j);
        
        psi(1) = (-4*k^3)/((exp(k*(T - tswitch1)) - 1)*((exp(k*(T - tswitch1)) - 1) * (1 - exp(-2*k*tswitch1)) - 2 * (exp(-k*tswitch1) - 1)));

        if psi(1) >= 0
            continue;
        end
        
        psi(2) = (psi(1)/k) * (exp(k*(T - tswitch1)) - 1);
        
        x21 =(1/(4*k))*(psi(2) - (psi(1)/k))*(exp(k*(t-tswitch1)) - exp(-k*(tswitch1 + t))) - (psi(1)/(2*k^2)) * (exp(-k * t) - 1);
        x2t1 =(1/(4*k))*(psi(2) - (psi(1)/k))*(1 - exp(-2*k*tswitch1)) - (psi(1)/(2*k^2)) * (exp(-k * tswitch1) - 1);
        x22 = x2t1*exp(k*(t - tswitch1)) + (psi(1)/(4*k^2)) * (exp(k*(T+t- 2*tswitch1)) - exp(k*(T-t)) + 1 - exp(k*(t - tswitch1)));
        
        
        psi21 = (psi(1)/k)*(exp(k*(T-t)) - 1);
        psi22= psi(2) * exp(-k*(t - tswitch1)) + (psi(1)/k) * (exp(k*(tswitch1 - t)) - 1);
        
        x21(t>=tswitch1) = x22(t>=tswitch1);
        psi21(t>=tswitch1) = psi22(t>=tswitch1);
        
        tVec1 = t(t<tswitch1);
        
        beforeSw = size(tVec1,2);
        x11 = zeros(size(tVec1));
        for i = 2:beforeSw
            x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
        end
        x11(1) = x11(2);
        
        tVec2 = t(t>=tswitch1);
        
        x12 = zeros(size(tVec2));
        for i = 2:size(tVec2,2)
            x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
        end
        x12(1) = x11(end);
        
        x1 = [x11, x12];
        
        u2 = zeros(1, NGrid);
        u2(t<tswitch1) = -k;
        u2(t>=tswitch1) = k;
        
        x1Matr(:,j + (l-1)*switchGrid) = x1;
        x2Matr(:,j+ (l-1)*switchGrid) = x21;
        psi2Matr(:,j+ (l-1)*switchGrid) = psi21;
        psi1Matr(:, j+ (l-1)*switchGrid) = psi(1);
        u2Matr(:,j + (l-1)*switchGrid) = u2;
    end
end


u1Matr = psi2Matr / 2;
IndVec = find(abs(x1Matr(end,:,1)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);

J9 = NaN;

if ~isempty(IndVec)
    
[J9, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 9;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:, minN);
u2MinMatr(:,r) = u2Matr(:, minN);
tSwitches = repmat(tSwitches, 1, 4);
SwitchTime(:,r) = [tSwitches(minN), NaN];
end

end

%% special 2

t = linspace(0,T,NGrid);
syms t1;
        x2T = eps;
        psi(2) = -2*eps^2;
        psi(1) = x2T * psi(2);
        tswitch1 = solve((1/(4*k))*(psi(2) - psi(1)/k) * (1 - exp(-2*k*t1)) - (psi(1)/(2*k^2)) * (exp(-k*t1) - 1) - eps);
        tswitch1 = double(tswitch1);
        ts1 = tswitch1(tswitch1>=0 & tswitch1<T);
        if ~isempty(ts1)
            ts1 = ts1(isreal(ts1));
        end
switchAm = size(ts1,2);
x1Matr = zeros(NGrid, 2*switchAm);
x2Matr = zeros(NGrid, 2*switchAm);
psi1Matr = zeros(NGrid, 2*switchAm);
psi2Matr = zeros(NGrid, 2*switchAm);

J10 = NaN;

edgesMatr =...
    [-M;M];
x2T = -eps;
psi1Vec = zeros(1,4);
if ~isempty(ts1)
for l = 1:2
    
    x10 = edgesMatr(l);
    
        for j = 1:size(t1,2)
            tswitch1 = ts1(j);
        x21 =(1/(4*k))*(psi(2) - (psi(1)/k))*(exp(k*(t-tswitch1)) - exp(-k*(tswitch1 + t))) - (psi(1)/(2*k^2)) * (exp(-k * t) - 1);
        psi21 = psi(2) * exp(k*(t - tswitch1)) + (psi(1)/k) * (1 - exp(-k*(tswitch1 - t)));
        
        x21(t>=tswitch1) = x2T;
        psi21(t>=tswitch1) = psi(2);
        
        
        tVec1 = t(t<tswitch1);
        
        beforeSw = size(tVec1,2);
        x11 = zeros(size(tVec1));
        for i = 2:beforeSw
            x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
        end
        if (size(x11,2) ~= 1)
            x11(1) = x11(2);
        else
            x11(1) = x12(1);
        end
        tVec2 = t(t>=tswitch1);
        
        x12 = zeros(size(tVec2));
        for i = 2:size(tVec2,2)
            x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
        end
        x12(1) = x11(end);
        
        x1 = [x11, x12];
        
        x1Matr(:,j + (l-1)) = x1;
        x2Matr(:,j+ (l-1)) = x21;
        psi2Matr(:,j+ (l - 1)) = psi21;
        psi1Matr(:, j+ (l - 1)) = psi(1);
        
        u2 = zeros(1, NGrid);
        u2(t<tswitch1) = -k;
        u2(t>=tswitch1) = k;

        u2Matr(:,j + (l-1)) = u2;
        end
end

u1Matr = psi2Matr / 2;
IndVec = find(abs(x1Matr(end,:)-L) < deltaEndSet & abs(x2Matr(end,:))<=eps);

J10 = NaN;

if ~isempty(IndVec)
    
[J10, minN] = min(trapz(t(1,:), u1Matr(:, IndVec) .^ 2));
minN = IndVec(minN);
r = 10;
if ~isempty(minN)
x1MinMatr(:, r) = x1Matr(:, minN);  
x2MinMatr(:, r) = x2Matr(:, minN);
psi1MinMatr(:, r) = psi1Matr(:, minN);
psi2MinMatr(:, r) = psi2Matr(:, minN);
u1MinMatr(:,r) = u1Matr(:, minN);
u2MinMatr(:,r) = u2Matr(:, minN);
tSwitches = repmat(tSwitches, 1, 4);
SwitchTime(:,r) = [tSwitches(minN), NaN];
end

end

end

%% The result N1

J = [J1, J2, J3, J4, J5];
FoundJInd = ~isnan(J);


if any(FoundJInd)
   [minJ, minN] = min(J(FoundJInd));
   r= find(FoundJInd);
   r = r(minN);
   
   figure;
   hold on;
   
   plot([-M, M],[0,0],'Color',[.35, .7, .2], 'LineWidth', 4);
   plot([L, L],[-eps,eps],'Color',[.7, .2, .35], 'LineWidth', 4);
   plot (x1MinMatr(:,r),x2MinMatr(:,r));
    plot(x1MinMatr(1,r), x2MinMatr(1,r), 'g.', 'MarkerSize', 20);
   plot(x1MinMatr(end,r), x2MinMatr(end,r), 'r.', 'MarkerSize', 20);
   xlabel('x_1');
   ylabel('x_2');
   
   hold off;
   figure;
   plot(t, psi1MinMatr(:,r));
   xlabel('t');
   ylabel('\psi_1');
   figure;
   plot(t, psi2MinMatr(:,r));
   xlabel('t');
   ylabel('\psi_2');
   figure;
   plot(t,u1MinMatr(:,r));
   xlabel('t');
   ylabel('u_1');
   figure;
   plot(t,u2MinMatr(:,r));
   xlabel('t');
   ylabel('u_2');
   disp('Switch times are (NaN means there was not a switch):');
   disp(SwitchTime(:,r));
else
    disp('---No  optimal solution found!!!---');
end

%% The result N2
J = [NaN, NaN, NaN, NaN, NaN, J6, J7, J8, J9, J10];

FoundJInd = ~isnan(J);
if any(FoundJInd)
   [minJ, minN] = min(J(FoundJInd));
   r= find(FoundJInd);
   r = r(minN);
   figure;
   hold on;
   
   plot([-M, M],[0,0],'Color',[.35, .7, .2], 'LineWidth', 4);
   plot([L, L],[-eps,eps],'Color',[.7, .2, .35], 'LineWidth', 4);
   plot (x1MinMatr(:,r),x2MinMatr(:,r));
    plot(x1MinMatr(1,r), x2MinMatr(1,r), 'g.', 'MarkerSize', 20);
   plot(x1MinMatr(end,r), x2MinMatr(end,r), 'r.', 'MarkerSize', 20);
   xlabel('x_1');
   ylabel('x_2');
   
   hold off;
   figure;
   plot(t, psi1MinMatr(:,r));
   xlabel('t');
   ylabel('\psi_1');
   figure;
   plot(t, psi2MinMatr(:,r));
   xlabel('t');
   ylabel('\psi_2');
   figure;
   plot(t,u1MinMatr(:,r));
   xlabel('t');
   ylabel('u_1');
   figure;
   plot(t,u2MinMatr(:,r));
   xlabel('t');
   ylabel('u_2');
   disp('Switch times are (NaN means there was not a switch):');
   disp(SwitchTime(:,r));
else
    disp('---No  optimal solution found!!!---');
end


function res = psi2x2Sys1(psi, tswitch1, x2T, T,k)
x2t1 =(1/(4*k))*(psi(2) + (psi(1)/k))*(exp(2*k*tswitch1) - 1) + (psi(1)/(2*k^2)) * (1 - exp(k * tswitch1));

res(1) = psi(1) + psi(2) * x2t1;
res(2) = -x2t1 + x2T*exp(-k*(tswitch1 - T)) + (1/(4*k)) * (psi(2) - (psi(1)/k)) * (1 - exp(2*k*(T - tswitch1))) - (psi(1)/(2*k^2)) * (exp(-k*(tswitch1 - T)) - 1);
end

function res = psi2x2Sys2(psi, tswitch1, x2T, T, k)
x2t1 =(1/(4*k))*(psi(2) - (psi(1)/k))*(1 - exp(-2*k*tswitch1)) - (psi(1)/(2*k^2)) * (exp(-k * tswitch1) - 1);

res(1) = psi(1) + psi(2) * x2t1;
res(2) = -x2t1 + x2T*exp(k*(tswitch1 - T)) + (1/(4*k)) * (psi(2) + (psi(1)/k)) * (exp(2*k*(tswitch1 - T)) - 1) + (psi(1)/(2*k^2)) * (1 - exp(k*(tswitch1 - T)));
end
