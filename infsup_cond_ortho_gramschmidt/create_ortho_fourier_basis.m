function [y,x] = create_ortho_fourier_basis(N,n)

x = linspace(0,1,n)';


y = x.^0;
plot(x,y)
hold on

for i = 1:N
    
    new_b = sin(x * pi * i);
    nbasis = size(y,2);
    
    for j = 1:nbasis
       new_b = new_b - trapz(x,y(:,j).*new_b)/trapz(x,y(:,j).*y(:,j)) * y(:,j);
    end
    new_b = new_b / sqrt(trapz(x,new_b.*new_b));
    
    y = [y new_b];
    plot(x,new_b)

    new_b = cos(x * pi * i);
    nbasis = size(y,2);
    
    for j = 1:nbasis
       new_b = new_b - trapz(x,y(:,j).*new_b)/trapz(x,y(:,j).*y(:,j)) * y(:,j);
    end
    new_b = new_b / sqrt(trapz(x,new_b.*new_b));
    
    y = [y new_b];
    %plot(x,new_b)
end

% nbasis = size(y,2);
% 
% M = zeros(nbasis);
% 
% for i = 1:nbasis
%     for j = 1:nbasis
%        M(i,j) = trapz(x,y(:,i).*y(:,j)); 
%     end
% end
