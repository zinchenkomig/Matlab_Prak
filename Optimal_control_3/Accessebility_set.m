alpha = 0.3;
t1 = 5;
t2 = 7;
N = 50; %frames for video
filename = 'one.avi';
clf;
reachsetdyn(alpha,t1,t2,N,filename);


function [x,y] = reachset(alpha, T)
eps = 0.01;
isPlotTraj = 0;
isPlotChange = 0;
isPlotStac = 0;
N = 20;
color = [0.1,0.4,0.6];
colorstart = [0.2, 0.8, 0.2];

options = odeset('Events', @SisZeroEvent);
[t1,x1,~,~,~] = ode45(@(t,x) s_small(t,x,alpha), [0,T], [0; 0], options);
[t2,x2,~,~,~] = ode45(@(t,x) s_small(t,x,-alpha), [0,T], [0; 0], options);
hold on;
axis equal;

if isPlotStac == 1
    syms x;
    sol = vpasolve(2*x^3 + sin(2*x) == alpha, x);
    x_stat1 = real(sol);
    plot(x_stat1, 0, 'Marker', '*', 'LineWidth', 1.5, 'Color', 'k');
    
    sol = vpasolve(2*x^3 + sin(2*x) == -alpha, x);
    x_stat2 = real(sol);
    plot(x_stat2,0, 'Marker', '*', 'LineWidth', 1.5, 'Color', 'k');
end


tt1 = linspace(0, t1(end), N);
tt2 = linspace(0, t2(end), N);
x11 = interp1(t1, x1(:,1), tt1);
x12 = interp1(t1, x1(:,2), tt1);
x21 = interp1(t2, x2(:,1), tt2);
x22 = interp1(t2, x2(:,2), tt2);

options = odeset('Events', @SisZeroEvent);

tEnd1 = zeros(1,N-1);
tEnd2 = zeros(1,N-1);

xEnd11 = zeros(1,N-1);
xEnd12 = zeros(1,N-1);
xEnd21 = zeros(1,N-1);
xEnd22 = zeros(1,N-1);

psiEnd1 = zeros(1,N-1);
psiEnd2 = zeros(1,N-1);

if isPlotChange == 1
    plot(x11,x12, 'r*');
    plot(x21,x22, 'r*');
end

for i = 1:N-1
    [t,xVec1,~,~,~] = ode45(@(t,x) S_Sys(t,x,-alpha), [tt1(i), T], [1;0;x11(i); x12(i)], options);
    if isPlotTraj == 1
        plot(xVec1(:,3), xVec1(:,4), 'Color',color);
    end
    tEnd1(i) = t(end);
    psiEnd1(i) = xVec1(end, 1);
    xEnd11(i) = xVec1(end,3);
    xEnd12(i) = xVec1(end,4);
    [t,xVec2,~,~,~] = ode45(@(t,x) S_Sys(t,x,alpha), [tt2(i), T], [1;0;x21(i); x22(i)], options);
    if isPlotTraj == 1
        plot(xVec2(:,3), xVec2(:,4), 'Color',color);
    end
    tEnd2(i) = t(end);
    psiEnd2(i) = xVec2(end, 1);
    xEnd21(i) = xVec2(end,3);
    xEnd22(i) = xVec2(end,4);
end

if isPlotChange == 1
    plot(xEnd11,xEnd12, 'r*');
    plot(xEnd21,xEnd22, 'r*');
end

if isPlotTraj == 1
    plot(x1(:,1), x1(:,2), 'Color',colorstart, 'LineWidth', 2);
    plot(x2(:,1), x2(:,2), 'Color',colorstart, 'LineWidth', 2);
end

while(any(tEnd1 < T-eps) || any(tEnd2 < T-eps))
    for i = 1:N-1
        if tEnd1(i) < T
            [t,xVec1,~,~,~] = ode45(@(t,x) S_Sys(t,x,alpha), [tEnd1(i), T], [psiEnd1(i);0;xEnd11(i); xEnd12(i)], options);
            if isPlotTraj == 1
                plot(xVec1(:,3), xVec1(:,4), 'Color',color);
            end
            tEnd1(i) = t(end);
            xEnd11(i) = xVec1(end,3);
            xEnd12(i) = xVec1(end,4);
            psiEnd1(i) = xVec1(end, 1);
        end
        if isPlotChange == 1
            plot(xEnd11,xEnd12, 'r*');
        end
        if tEnd2(i) < T
            [t,xVec2,~,~,~] = ode45(@(t,x) S_Sys(t,x,-alpha), [tEnd2(i), T], [psiEnd2(i);0;xEnd21(i); xEnd22(i)], options);
            if isPlotTraj == 1
                plot(xVec2(:,3), xVec2(:,4), 'Color',color);
            end
            tEnd2(i) = t(end);
            xEnd21(i) = xVec2(end,3);
            xEnd22(i) = xVec2(end,4);
            psiEnd2(i) = xVec2(end, 1);
        end
        if isPlotChange == 1
            plot(xEnd21,xEnd22, 'r*');
        end
    end
    alpha = - alpha;
end

x = [xEnd11,xEnd21,xEnd11(1)];
y = [xEnd12,xEnd22,xEnd12(1)];

    function res = s_small(t,x, alpha)
        res = [x(2); -2*x(1)^3 * cos(x(2)) - sin(2*x(1)) - x(2) + alpha];
    end

    function res = S_Sys(t,y, alpha)
        res = [-y(2) * 6 * y(3) ^ 2 * cos(y(4)) - 2 * y(2) * cos(2 * y(3));...
            y(1) + y(2) * 2 * y(3)^3 * sin(y(4)) - y(2);...
            y(4);...
            -2*y(3)^3 * cos(y(4)) - sin(2*y(3)) - y(4) + alpha];
    end

    function [position, isterminal, direction] = SisZeroEvent(t,x)
        position = x(2);
        isterminal = 1;
        direction = 0;
    end
end


function reachsetdyn(alpha, t1, t2, N, filename)
colorend = [0.8, 0.25, 0.25];
if (t1 ~= t2)
    TVec = linspace (t1, t2, N);
    v = VideoWriter(filename, 'Motion JPEG AVI');
    v.FrameRate = 10;
    open(v);
    for i = 1:N
        clf;
        [x,y] = reachset(alpha, TVec(i));
        plot(x,y,'Color',colorend, 'LineWidth', 2);
        frames(i) = getframe(gcf);
        writeVideo(v,frames(i));
    end
    close(v);
else
   [x,y] = reachset(alpha, t1);
   plot(x,y, 'Color',colorend, 'LineWidth', 2);
end
end