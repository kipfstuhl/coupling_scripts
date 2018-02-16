function [F,J] = evaluate_function_and_jacobian(u)
    global mat jac rhs
    F = mat(u)*u-rhs;
    if nargout > 1
        J = jac(u);
    end
end