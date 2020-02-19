% test steepest descent method

% 1 Dim
r = rand;
f = @(x) r*(23*x.^2+ x^5 - 4);

df = @(x)r*(2*23*x+ 5*x^4) ;

a = gradient_descent(f, df, [1]);
f(a); % works damn well

% 2 Dim
f = @(x) [x(1)^2-x(2)
    x(2)-4];
df = @(x) [ 2*x(1) , -1
    0, 1];

a = gradient_descent(f, df, rand*[1,2]);
f(a);

% 2 Dim - harder 
r = rand;
f = @(x) [x(1)^2-x(2)+x(2)*x(1)^5
    sin(x(2))+x(2)*x(1)^5];
df = @(x) [ 2*x(1)+x(2)*5*x(1)^4 , -1+x(1)^5
    +x(2)*5*x(1)^4, cos(x(2))+x(1)^5];
gradf = @(x) 2*f(x)'*df(x);

% this works
a = gradient_descent(f, df, [1,2])
f(a)
