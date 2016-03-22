function [ F ] = state_function( alpha, beta, u )
%STATE_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

len = length(alpha);
if length(alpha)~=length(beta)
    error('The length of alpha parameter doesn''t match the length of beta parameter');
end;

F=zeros(len,len);
step_a= alpha(1,2)- alpha(1,1);
step_b= beta(2,1)- beta(1,1);

for i=1:1:len
    for j=i:1:len
        for k=i:1:j
            for l=j:-1:i
                F(i,j)=F(i,j)+u(k,l);
            end;
        end;
    end;
end;

F=F.*(step_a*step_b);
end