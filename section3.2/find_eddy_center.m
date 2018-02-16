function xeddy = find_eddy_center(sol,x0)

fun = @(x) evaluate_solution_velocity_norm(sol,x);

xeddy = fminunc(fun,x0);

end

function res = evaluate_solution_velocity_norm(sol,x)
    res = sqrt(evaluate_fe_function(sol.u1,sol.fespace_u,x)^2 + evaluate_fe_function(sol.u2,sol.fespace_u,x)^2);

    if (res == 0)
        res =  Inf;
    end
end
