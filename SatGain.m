function g = sat_gain(Ew,w,g0,Is,Wp)
Imean = mean(abs(Ew).^2); % mean roundtrip intensity
Wmin = min(w);
if Wmin < Wp
  Iss = @(waistMin) Is*waistMin.^2/Wp.^2;
  g = g0*1/(1 + Imean/Iss(Wmin)); 
else
  g = g0*1/(1 + Imean/Is); 
end
end