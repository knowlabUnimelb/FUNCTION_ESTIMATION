function y = gety(x, p, noise, type)

% mu = 0; sd = noise;
% adj = normrnd(mu, sd, size(x));
%
% y = polyval(p, x) + adj; % Generate y values adjusted by normally distributed noise

switch type
    case 'poly'
        mu = 0; sd = noise;
        adj = normrnd(mu, sd, size(x));
        
        y = polyval(p(~isnan(p)), x);
%         y = y + (std(y) * adj); % Generate y values adjusted by normally distributed noise
        y = y + adj;
        
    case 'sim'
        n = .75;
        tx = linspace(-6, 6, 10);
        
        mu = 0; sd = n;
        adj = normrnd(mu, sd, size(tx));
        
        ty = polyval(p(~isnan(p)), tx);
        ty = ty + (std(ty) * adj); % Generate y data with noise
        
        % Fit 7th degree poly to data
        pco = polyfit(tx,ty, 7);
        
        % Now, generate real data
        mu = 0; sd = noise;
        adj = normrnd(mu, sd, size(x));
        
        y = polyval(pco(~isnan(pco)), x);
        y = y + (std(y) * adj); % Generat
end