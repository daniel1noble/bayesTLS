

l = 0.05
u = 0.98
k = 8
time_min = c(1, 5, 10, 20, 40, 60, 120, 180)
mid = log10(60) # Note that because log10 transformed time_min also needs to be log10 so on same units

pij = l + ((u-l) / (1 + exp(k * (log10(time_min) - mid))))


plot(pij ~ time_min, type = "l")