function [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,xi,eta)
if elem_type == 1
    D_N = 0.25*[eta-1 1-eta 1+eta -1-eta; xi-1 -xi-1 xi+1 1-xi];
elseif elem_type == 2
    D_N = [ - (xi/4 - 1/4)*(eta - 1) - ((eta - 1)*(eta + xi + 1))/4, ((eta - 1)*(eta - xi + 1))/4 - (xi/4 + 1/4)*(eta - 1), (xi/4 + 1/4)*(eta + 1) + ((eta + 1)*(eta + xi - 1))/4, (xi/4 - 1/4)*(eta + 1) + ((eta + 1)*(xi - eta + 1))/4, (eta/2 - 1/2)*(xi - 1) + (eta/2 - 1/2)*(xi + 1),                          -((eta - 1)*(eta + 1))/2, - (xi/2 + 1/2)*(eta + 1) - ((eta + 1)*(xi - 1))/2,                         ((eta - 1)*(eta + 1))/2;...
- (xi/4 - 1/4)*(eta - 1) - (xi/4 - 1/4)*(eta + xi + 1),  (xi/4 + 1/4)*(eta - xi + 1) + (xi/4 + 1/4)*(eta - 1),  (xi/4 + 1/4)*(eta + 1) + (xi/4 + 1/4)*(eta + xi - 1),  (xi/4 - 1/4)*(xi - eta + 1) - (xi/4 - 1/4)*(eta + 1),                           ((xi - 1)*(xi + 1))/2, - (xi/2 + 1/2)*(eta - 1) - (xi/2 + 1/2)*(eta + 1),                            -(xi/2 + 1/2)*(xi - 1), (xi/2 - 1/2)*(eta - 1) + (xi/2 - 1/2)*(eta + 1)];
end
% size(D_N)
J = D_N*[X Y];
Jinv = inv(J);
jac = det(J);
B = J\D_N;
end