function [ alpha, beta, u ] = norm_ab( max_ab, step, mean_a, mean_b,sigma_a, sigma_b)
%NORM_AB creates a normal distribution of beta and alpha parameters of
%hysterons whose left border is alpha and right border is beta

mu = [mean_a mean_b];
sigma = [sigma_a 0; 0 sigma_b];
x1_for_alpha = -max_ab:step:max_ab; 
x2_for_beta = -max_ab:step:max_ab; 
[alpha, beta] = meshgrid(x1_for_alpha,x2_for_beta);
u = mvnpdf([alpha(:) beta(:)],mu,sigma);
u = reshape(u,length(x2_for_beta),length(x1_for_alpha));

for i=1:1:length(x1_for_alpha)
    for j=1:1:length(x2_for_beta)
        if alpha(i,j)<beta(i,j)
            u(i,j)=0;
        end;
    end;
end;
surf(x1_for_alpha, x2_for_beta,u);
caxis([min(u(:))-.5*range(u(:)),max(u(:))]);
axis([-max_ab max_ab -max_ab max_ab 0 .4])
xlabel('x1 for alpha'); ylabel('x2 for beta'); zlabel('Probability Density');

end

