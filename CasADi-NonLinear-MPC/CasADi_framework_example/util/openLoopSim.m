function [T, X] = openLoopSim( ffun, Times, x0, U, D, p )



N = length(Times)-1;

xk = x0;

T = Times(1);
X = x0';

% options.N = 20;

for i = 1:N
    
    [Tode, Xode] = ode45(  ffun, [Times(i) , Times(i+1)], xk, [], U(:,i), D(:,i), p );
    
    T = [T ; Tode(2:end)];
    X = [X ; Xode(2:end,:)];

    xk = X(end,:)';
    
end

end