function l = likelihood(x, y, x_j, y_j, t_obs, v, sigma)
l = 1;
for i = 1:length(x_j)
    t_j = 1/v * sqrt((x_j(i) - x).^2 + (y_j(i) - y).^2);
    l_aux = (1/(sqrt(2*pi)*sigma)) *exp(-((t_j - t_obs(i)).^2)/(2*sigma));
    l = l .* l_aux;
end
end