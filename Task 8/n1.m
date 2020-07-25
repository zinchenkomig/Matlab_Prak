%Monte-Carlo task 8
n = 40;
MinVisits = 400;

x = linspace(-1,1,n);
[XGrid, YGrid] = meshgrid(x);
IsInD = (XGrid.^2 + YGrid.^2 <= 1);

NeighboursCount = [IsInD(2:end, :); zeros(1, size(IsInD, 2))] +...
    [zeros(1, size(IsInD, 2)); IsInD(1:end-1, :)] +...
    [IsInD(:, 2:end), zeros(size(IsInD, 1), 1)] +...
    [zeros(size(IsInD, 1),1), IsInD(:, 1:end-1)];

NeighboursCount(~IsInD) = 0;
IsEdge = (NeighboursCount < 4 & NeighboursCount > 0);
Sums = zeros(size(XGrid));
Sums(IsEdge) = f(XGrid(IsEdge), YGrid(IsEdge));

VisitsAmount = zeros(size(XGrid));

[i , j] = find(IsEdge, 1);
S = Sums(i,j);

while min(VisitsAmount(IsInD & ~IsEdge)) < MinVisits
    r = randi(4);
    switch r
        case 1 %up
            if i>1
                i = i - 1;
            end
            
        case 2 %down
            if i<n
                i = i + 1;
            end
            
        case 3 %left
            
            if j>1
                j = j - 1;
            end
            
        case 4 %right
            
            if j<n
                j = j + 1;
            end
        otherwise
            disp('Error');
    end
    
    if IsInD(i, j)
        if ~IsEdge(i,j)
            Sums(i, j) = Sums(i, j) + S;
            VisitsAmount(i, j) = VisitsAmount(i,j) + 1;
        else
            S = Sums(i, j);
        end
    end
    
end

VisitsAmount(IsEdge) = 1;
Sums = Sums ./ VisitsAmount;

figure;
surf(XGrid, YGrid, Sums);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');

figure;
u = f(XGrid, YGrid);
u(~IsInD) = NaN;
surf(XGrid, YGrid, u);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');

function z = f(x,y)
z = x.^2 - y.^2;
end