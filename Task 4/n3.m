nVec = 1:1000:1e5;
averaging = 20;
k = size(nVec, 2);
times_pairs = zeros(1,k);
times_phone = zeros(1,k);
times_avg = zeros(1,k);
disp('Please wait...');
for i = 1:k
    n = nVec(i);
    for t = 1:averaging
        tstart = tic;
        a = standard_pairs(n);
        times_avg (t) = toc(tstart);
    end
    times_pairs(i) = mean(times_avg);
    
    for t = 1:averaging
        tstart = tic;
        b = standard_phoneNewMan(n);
        times_avg(t) = toc(tstart);
    end
    times_phone(i) = mean(times_avg);
end

plot(nVec, times_phone);
hold on;
plot(nVec, times_pairs);
legend('fonNewMan', 'Pairs');
xlabel('sample size');
ylabel('time');
