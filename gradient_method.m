function [x,fval]=gradient_method(func,x0,c,epsilon)
% Gradient method with backtracking stepsize rule
%
% INPUT
%=======================================
% func .......... must output objective function value and gradient of the objective function
% x0......... initial point
% c ..... tolerance parameter for the sufficient descent (Armijo)
% condition
% epsilon ... tolerance parameter for stopping rule
% OUTPUT
%=======================================
% x ......... optimal solution (up to a tolerance) 
%             of min f(x)
% fun_val ... optimal function value
beta = 0.5; %parameter to reduce alpha; could be some other value between 0 and 1
x=x0;
[fval, grad]=func(x);
iter=0;
fprintf('iter= %2d x = %s \tf(x)=%10.10f\n',iter,sprintf('%1.4f ',x'),fval) 
plot(x(1),x(2),'r*','MarkerSize',16,'LineWidth',2);
while (norm(grad)>epsilon)
    iter=iter+1;
    alpha=1;
    d = -grad; %search direction
    while true        
        x_new = x + alpha*d;
        if (func(x_new) - fval <= c*alpha*d'*grad) %if sufficient decrease is satisfied, exit
            break;
        end
        alpha=beta*alpha; %otherwise, reduce step size
    end
    x=x_new;
    plot(x(1),x(2),'b*','MarkerSize',16,'LineWidth',2);
   [fval, grad]=func(x);
fprintf('iter= %2d x = %s \tf(x)=%10.10f\n',iter,sprintf('%1.4f ',x'),fval) 
end
