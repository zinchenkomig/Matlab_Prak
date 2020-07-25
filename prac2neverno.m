% %psi1 < 0; psi2 <= 0
% %
% T = 5;
% k = 1;
% M = 24;
% L = 13;
% eps = 7;
% NGrid = 1000;
% BeginN = 2;
% switchGrid = 20;
% deltaEndSet = 0.3;




% tVec = linspace(0, T, switchGrid);
% tswitch1 = tVec(2:end);
% tVec = tswitch1;
%
% A = repelem(tVec, 1, size(tVec, 2))';
% B = repmat(tVec, 1, size(tVec, 2))';
% C = A < B;
% tswitch1 = A(C);
% tswitch2 = B(C);
%
% swSize = size(tswitch1,1);
%
% psi1 = 4*k^3 ./ ((2 - exp(k*(tswitch1 - tswitch2)) - exp(k*(tswitch2 - tswitch1))) .* (1 - exp(k*(tswitch2 - tswitch1))));
% %all possible psi1 for switch moments
%
%
% t = linspace(0,T,NGrid)';
% t = repmat(t,1,size(tswitch1,1))';
%
% psi21 = (psi1/k).*(1 - exp(k*(t - tswitch1)));
%
%
% psi2t2 = (psi1/k) .* (1 - exp(k*(tswitch2 - tswitch1)));
% psi22 = exp(k*(tswitch2 - t)) .* (psi2t2 + (psi1/k)) - psi1/k;
% psi21(t>=tswitch2) = psi22 (t>=tswitch2);
%
%
% plot(t(1,:)',psi21(1,:)');
%
%
% x21 = zeros(size(t));
% x22 = (psi1/(4*k^2)) .* (2 - exp(k*(tswitch1-t)) - exp(k*(t - tswitch1)));
% x2t2 = (psi1/(4*k^2)) .* (2 - exp(k*(tswitch1-tswitch2)) - exp(k*(tswitch2 - tswitch1)));
%
% x23 = x2t2 .* exp(k*(t - tswitch2)) + (1/(2*k)) * sinh(k*(t-tswitch2)) .* (psi2t2 + psi1/k) +...
%     (psi1/(2*k^2)) .* (1 - exp(k*(t - tswitch2)));
%
% x21(t>tswitch1 & t<tswitch2) = x22(t>tswitch1 & t<tswitch2);
% x21(t>=tswitch2) = x23(t>=tswitch2);
% %plot(t',x21');
% %plot(t(93,:), x21(93,:));
%
% x10 = linspace(-M, M, BeginN);
% x11 = zeros(swSize, NGrid, BeginN);
%
% x1t1 = zeros(swSize,2);
% x1t2 = zeros(swSize,2);
% for i = 1:BeginN
%     xbuf1 = x10(i) - k*t;
%     x1t1(:,i) = x10(i)-k*tswitch1;
%
%     xbuf2 = x1t1(:,i) + ...
%         (psi1/(4*k^2)) .* (2*(t-tswitch1) + (2/k) * sinh(k*(tswitch1- t))) - k*(t-tswitch1);
%
%     x1t2(:,i) = x1t1(:,i) + ...
%         (psi1/(4*k^2)) .* (2*(tswitch2-tswitch1) + (2/k) * sinh(k*(tswitch1- tswitch2))) - k*(tswitch2-tswitch1);
%
%     xbuf3 = x1t2(:,i) + ...
%        (x2t2/k).* (exp(k*(t - tswitch2)) - 1) + ...
%        (1/(2*k^2)) * (psi2t2 + psi1/k).*(cosh(k*(t-tswitch2))  - 1) +...
%        (psi1/(2 * k^2)) .* (t-tswitch2) -...
%        (psi1/(2*k^3)) .* (exp(k*(t - tswitch1)) - exp(k*(tswitch2 - tswitch1)))+...
%        k * (t - tswitch2);
%
%    xbuf1(t>tswitch1 & t<tswitch2) = xbuf2(t>tswitch1 & t<tswitch2);
%    xbuf1(t>=tswitch2) = xbuf3(t>=tswitch2);
%    x11(:,:,i) = xbuf1;
% end
%
% IndVec1 = find(abs(x11(:,end,1)-L) < deltaEndSet & abs(x21(:,end))<eps);
% IndVec2 = find(abs(x11(:,end,2)-L) < deltaEndSet & abs(x21(:,end))<eps);
%
% TrajN = 66;
% xPlMat = x11(:, :, 2)';
% yPlMat = x21(:, :)';
%
% % xPlMat = [xPlMat, x11(IndVec2, :, 2)'];
% % yPlMat = [yPlMat, x21(IndVec2, :)'];
% % IndVec = [IndVec1, IndVec2];
% %plot(xPlMat(:,:), yPlMat(:,:));
%
% psi21(t>tswitch1) = 0;
% psi21 = psi21';
%
% %J = trapz(t(1, :),(psi21(:, IndVec)).^2 /2);



% %Transvers: psi2(T) = 0; x2(T) != +-eps;
% edgesMatr = [-M; M];
% 
% for l = 1:2
%     
%     x10 = edgesMatr(l);
%     
%     for j = 2:switchGrid
%         
%         tswitch1 = tSwitches(j);
%         
%         psi(1) = (-2*k^2)/((1 - exp(k*(tswitch1 - T)))*(tswitch1 * (2 - exp(k*(tswitch1 - T))) + (1/k) * (exp(-k * tswitch1) - 1)));
% 
%         if psi(1) <= 0
%             break;
%         end
%         
%         psi(2) = (psi(1)/k) * (1 - exp(k*(t - T)));
%         
%         x21 = (t/2) .* exp(k*(tswitch1 - t))*(psi(2) + (psi(1)/k)) + (psi(1)/(2*k^2)) * (exp(-k*t) - 1);
%         x2t1 = (tswitch1/2) * (psi(2) + (psi(1)/k)) + (psi(1)/(2*k^2)) * (exp(-k*tswitch1) - 1);        
%         x22 = x2t1 * exp(-k*(t - tswitch1)) + (psi(1)/(4*k^2)) * (exp(k*(2*tswitch1 - T - t)) - exp(k*(t-T)) + 2 - 2 * exp(k*(tswitch1 - t)));       
%         
%         psi21 = psi(2) * exp(-k*(t - tswitch1)) + (psi(1)/k) * (exp(k*(tswitch1 - t)) - 1);
%         psi22 = (psi(1)/k) * (1 - exp(k*(t - T)));
%         
%         x21(t>=tswitch1) = x22(t>=tswitch1);
%         psi21(t>=tswitch1) = psi22(t>=tswitch1);
%         
%         
%         tVec1 = t(t<tswitch1);
%         
%         beforeSw = size(tVec1,2);
%         x11 = zeros(size(tVec1));
%         for i = 2:beforeSw
%             x11(i) = x10 + trapz(tVec1(1:i), x21(1:i) - k);
%         end
%         x11(1) = x11(2);
%         
%         tVec2 = t(t>=tswitch1);
%         
%         x12 = zeros(size(tVec2));
%         for i = 2:size(tVec2,2)
%             x12(i) = x11(end) + trapz(tVec2(1:i), x21(beforeSw+1:beforeSw + i) + k);
%         end
%         x12(1) = x11(end);
%         
%         x1 = [x11, x12];
%         
%         x1Matr(:,j + l*switchGrid) = x1;
%         x2Matr(:,j+ l*switchGrid) = x21;
%         psi2Matr(:,j+ l*switchGrid) = psi21;
%         psi1Matr(:, j+ l*switchGrid) = psi(1);
%     end
% end
% 
% plot(x1Matr, x2Matr);

%% 3
% psi1<0; psi2<=0;

% One switch

% T = 8;
% k = 1;
% M = 11;
% L = 16.93;
% eps = 20;
% NGrid = 1000;
% BeginN = 2;
% switchGrid = 20;
% deltaEndSet = 0.3;

tVec = linspace(0, T, switchGrid);
tswitch1 = tVec(2:end)';

psi1 = 4* eps * k^2 ./ (2 - exp(k*(tswitch1 - T)) - exp(k*(T - tswitch1)));

%all possible psi1 for switch moments
swSize = size(tswitch1,1);

t = linspace(0,T,NGrid)';
t = repmat(t,1,swSize)';

psi21 = (psi1/k) .* (1 - exp(k*(t - tswitch1)));


x21 = zeros(size(t));
x22 = (psi1/(4*k^2)) .* (2 - exp(k*(tswitch1 - t)) - exp(k*(t - tswitch1)));
x21(t>=tswitch1) = x22(t>=tswitch1);
x21 = x21';



x10 = linspace(-M, M, BeginN);
x11 = zeros(NGrid, swSize, BeginN);
x1t1 = zeros(swSize,2);
x1t2 = zeros(swSize,2);

for i = 1:BeginN
    xbuf1 = x10(i) - k*t;
    x1t1(:,i) = x10(i) - k*tswitch1;
    xbuf2 = x1t1(:,i) + ...
        (psi1/(4*k^2)) .* (2*(t-tswitch1) + (2/k) * sinh(k*(tswitch1- t))) - k*(t-tswitch1);
    
    xbuf1(t>=tswitch1) = xbuf2(t>=tswitch1);
    x11(:,:,i) = xbuf1';
    
end

IndVec1 = find(abs(x11(end,:,1)-L) < deltaEndSet & abs(x21(end,:))<=eps);
IndVec2 = find(abs(x11(end,:,2)-L) < deltaEndSet & abs(x21(end,:))<=eps);

psi21 = psi21';

u1Matr = zeros(switchGrid,NGrid);
u1Matr(t>=tswitch1) = psi21(t>=tswitch1) / 2;

x1Matr = x11(:, :, 1);
x2Matr = x21(:, :);

psi2Matr = psi21;
psi1Matr = repmat(psi1', NGrid, 1);

x1Matr = [x1Matr, x11(:, :, 2)];
x2Matr = [x2Matr, x21(:, :)];
IndVec = [IndVec1, IndVec2];


J3 = NaN;
if ~isempty(IndVec)
[J3, minN] = min(trapz(t(1,:), u1Matr(:, IndVec)' .^ 2));
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

