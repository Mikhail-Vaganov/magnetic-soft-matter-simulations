function [ alpha, beta, u ] = random_ab( max_ab, step)
%RANDOM_AB creates a random distribution of beta and alpha parameters of
%hysterons whose left border is alpha and right border is beta

x1_for_alpha = -max_ab:step:max_ab; 
x2_for_beta = -max_ab:step:max_ab; 
[alpha, beta] = meshgrid(x1_for_alpha,x2_for_beta);

for i=1:1:length(x1_for_alpha)
    for j=1:1:length(x2_for_beta)
        if alpha(i,j)<beta(i,j)
            u(i,j)=0;
        else
            u(i,j)=rand(1);
        end;
    end;
end;

surf(x1_for_alpha, x2_for_beta,u);
caxis([min(u(:))-.5*range(u(:)),max(u(:))]);
axis([-max_ab max_ab -max_ab max_ab 0 .4])
xlabel('x1 for alpha'); ylabel('x2 for beta'); zlabel('Probability Density');

end

