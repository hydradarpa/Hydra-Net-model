

figure;
plot((1:t)*dt,I_repr_ext,(1:t)*dt,I_repr_EE);
set(gca,'FontSize',20);             % set the axis with big font
title('Synaptic input to the representative neuron');
box('off');
xlabel('time (ms)');
ylabel('Input (\muA)');