function X = TTorthogonalizeLR(X)
% TTORTHOGONALIZELR  left to right orthogonalization of a TT-tensor.
%
%  X = TTORTHOGONALIZELR(Y) returns a TT-tensor X  with left-orthogonal cores X{1}, ..., X{N-1}.

[N,I,~] = TTsizes(X);

for n = 1 : N - 1
    [X{n}, R] = qr(X{n}, 0);
    X{n + 1} = h2v(R * v2h(X{n + 1}, I(n + 1)), I(n + 1));
end

