function [ ACOnextroute ] = rouletteWheel( P )

cumsumP = cumsum(P);

r = rand();

ACOnextroute = find(r <= cumsumP);

ACOnextroute = ACOnextroute(1);

end