function k = makePolyKernel(x, xprime, sigma, p1, p2, simple)

switch simple
    case 1
        for i = 1:length(x)
            for j = 1:length(xprime)
                k(i,j) = gensimplepoly(x(i), p1)' *  sigma * gensimplepoly(xprime(j), p2);
            end
        end
    case 2
        for i = 1:length(x)
            for j = 1:length(xprime)
                k(i,j) =  genorthopoly(x(i), p1)' *  sigma * genorthopoly(xprime(j), p2);
            end
        end
end