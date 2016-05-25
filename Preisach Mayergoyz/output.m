function [ output ] = output(input, t, F, alpha, beta)
%OUTPUT Summary of this function goes here
%   Detailed explanation goes here

apprF = fit([alpha(:),beta(:)],F(:),'poly23');

len = size(alpha,2);

max_alpha = max(max(alpha));
min_beta = min(min(beta));

for i=1:1:length(t)
    
    if(input(i)>max_alpha)
       output(i) = apprF(max_alpha,min_beta);
       continue;
    end;
    if(input(i)<min_beta)
       output(i) = -apprF(max_alpha,min_beta);
       continue;
    end;
    
    output(i) = -apprF(max_alpha,min_beta);
    input_t(i)=input(i);
    
    
    
    
    maximum = max_alpha;
    minimum = min_beta;
    t_maximum=t(1);
    t_minimum=t(1);
    if i>2        
        [max_values, t_max] = findpeaks(input_t, t(1:i));
        [min_values, t_min] = findpeaks(-input_t,t(1:i));
        maximum=[maximum,max_values];
        minimum=[minimum,min_values];
        t_maximum=[t_maximum,t_max];
        t_minimum=[t_minimum,t_min];
    end
        
    n=min(length(maximum), length(minimum));
        
    for k=2:1:n
        output(i)=output(i)+2*(apprF(maximum(k), minimum(k-1)-apprF(maximum(k),minimum(k))));
    end
    
    if i==1
        output(i)=output(i)+2*(apprF(input(i), minimum(1)));
    else if input(i) < input(i-1)
        if(n==1)
            output(i)=output(i)+2*(apprF(maximum(n), minimum(1))-apprF(maximum(n),input(i)));
        else
            output(i)=output(i)+2*(apprF(maximum(n), minimum(n-1))-apprF(maximum(n),input(i)));
        end;
    else
        if(n==1)
            output(i)=output(i)+2*(apprF(input(i), minimum(1)));
        else
            output(i)=output(i)+2*(apprF(input(i), minimum(n-1)));
        end
    end
end;



end

