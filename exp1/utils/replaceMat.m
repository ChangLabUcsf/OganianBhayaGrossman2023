function [B] = replaceMat(A, orig, subs)
    B = A;
    B(B==orig) = subs;
end