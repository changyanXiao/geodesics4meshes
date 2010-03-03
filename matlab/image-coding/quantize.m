function [vI,vQ] = quantize(v,T)

% quantized integer values
vI = floor(abs(v/T)).*sign(v);
% de-quantized values from vI
vQ = sign(vI) .* (abs(vI)+.5) * T;