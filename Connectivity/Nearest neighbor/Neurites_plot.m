
%%
figure
Q=sum(A);
histogram(Q)
box off
set(gca,'Fontsize',30)
ylabel('Count')
xlabel('Number of connections')
title('Number of neurites in CB network')