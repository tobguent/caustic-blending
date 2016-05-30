function [ T ] = pmInterpolatePnts( A, I, B, t, w )
% computes a cubic chordal-parameterized b-spline going through the points A P and B, whereas P
% is computed from a weighting between T and the chordal blend between A and B.

delta1 = sqrt(sum((A-I).^2));
delta2 = sqrt(sum((I-B).^2));

r = delta1 ./ (delta1+delta2);
rr = repmat(r,size(A,1),1);

P = A + (B-A) .* rr;    % point between A and B
Q = P + (I-P) * w;      % point between the linear blend and I

T = pmInterpolatePnts_mex(A, Q, B, t, r);

end

