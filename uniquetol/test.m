foo = @(x) plus(1, x * 1e-14);

a = arrayfun(foo, -107:108);
x = a([1, 108, end]);
y = a([1, 100]);
z = a([100, end]);

uniquetol(a)
uniquetol(x)
uniquetol(y)
uniquetol(z)

% a0 = arrayfun(foo, -49:50);
% e0 = a0([1, 50, end]);
% uniquetol(a0)
% uniquetol(e0)

% a1 = arrayfun(foo, -99:100);
% e1 = a1([1, 100, end]);
% uniquetol(a1)
% uniquetol(e1)