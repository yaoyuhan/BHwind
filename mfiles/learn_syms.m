clear all

syms x y real
[x, y] = solve('x^2 + x*y + y = 3', 'x^2 - x*4 + 3 = 0');
solution  = [x, y]