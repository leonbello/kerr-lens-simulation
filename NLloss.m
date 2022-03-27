function l = loss_gain(w,Wp)
for i = 1:length(w)
    if w(i) <= Wp
        l(i) = 1;
    else
        loss_fun = @(w) 1./(1+(w-Wp).^2/(2*(Wp).^2));
        l(i) = loss_fun(w(i));
    end
end
end