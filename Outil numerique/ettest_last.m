%% étude de l'histogramme
function E = ettest_last(residu,paramDist,critere)
    warning('off')
    Snoise = residu;
    mu = paramDist(1);
    sigma = paramDist(2);
    
    nbins = 50;
    [nelements,bincenters] = hist(Snoise,nbins);
    X = bincenters;
    Y = nelements;
    lengthBin = X(3)-X(2);
    %scalingFactor = size(Snoise,1) * lengthBin;

    %xx = linspace(min(X)-lengthBin/2,max(X)+lengthBin/2,350);
    xx = linspace(mu-3*sigma,mu+3*sigma,350);
    yy = normpdf(xx,mu, sigma);

    yC = chebfun(yy',[xx(1) xx(end)],'equi');
    NN = length(Snoise);
    Nb = 0;
    hh=1;
    
    for j=1:nbins
        a = X(j)-lengthBin/2;
        b = X(j)+lengthBin/2;
        x = chebfun('x',[a b]);
        p(j) = sum(yC(x));

        mumu(j) = NN*p(j);
        sigsig(j) = sqrt(NN*p(j)*(1-p(j)));

        % pourcentage d'efficacité de l'histogramme (doit être 68%)
        if mumu(j)-sigsig(j)/hh<=Y(j) && Y(j)<=mumu(j)+sigsig(j)/hh
            Nb=Nb+1;
        end
    end
    
    E=Nb./nbins;
    
    %if E>=0.68
    if E>=critere
        R=1;
    else
        R=0;
    end
    
end