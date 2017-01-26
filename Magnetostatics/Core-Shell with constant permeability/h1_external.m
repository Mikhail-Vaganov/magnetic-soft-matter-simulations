function [ output] = h1_external( chi, h0, f, hm )

output = (2*(f^3-1)*chi*(hm*(3+chi)-3*h0*chi))/(3*(-2*chi^2+f^3*(9+9*chi+2*chi^2)))+h0;

end

