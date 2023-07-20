ns = 500:500:8000;
ts = 0*ns;
for j = 1:length(ns)
    A = randn(ns(j),ns(j)); 
    % create a matrix with random number entries
    x = randn(ns(j),1);
    tic; % start the timer
    for i = 1:15 % repeat fifteen times
        A*x;
    end
    t = toc; % read the timer
    ts(j) = t/15;
end

loglog(ns,ts,'.-r','MarkerSize',25,'LineWidth',1.5) % plot the data
xlabel('matrix size (n)')
ylabel('time (sec)')
hold on
set(gca,'FontSize',16);
loglog(ns,ts(end)*(ns/ns(end)).^2,'--k','LineWidth',1.5) 
% plot a reference line

grid on
axis tight
legend('data','O(nˆ2)','location','southeast')