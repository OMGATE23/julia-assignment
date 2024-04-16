# function trapezoidal(f::Function , limits::Tuple{Real , Real} , n::Int64)
#     a,b = limits
#     delta_x = (b-a)/n
#     xi = a:delta_x:b
#     return (delta_x/2)*(f(xi[1]) + f(xi[end] + 2 * sum(f.(xi[2:end]))))
# end

# function_question_2 = x -> sqrt(abs(x))log(x)

# trapezoidal(function_question_2 , (0, 1) , 10)

using Plots

function trapezoidal(f::Function, limits::Tuple{Real, Real}, n::Int64)
    a, b = limits
    delta_x = (b - a) / n
    xi = a:delta_x:b
    sum_val = sum(f.(xi[2:end-1])) + (f(a) + f(b)) / 2 
    return delta_x * sum_val
end

function function_question_2(x)
    if x == 0
        return 0  # handle edge case for 0 value
    else
        return sqrt(abs(x)) * log(x)
    end
end

true_value = -4/9

n_values = 2 .^ (1:10)
errors = []
for n in n_values
    integral_approx = trapezoidal(function_question_2, (0, 1), n)
    println(integral_approx)
    error = abs(integral_approx - true_value)
    push!(errors, error)
end
plot(n_values, errors, xscale=:log10, yscale=:log10, marker=:circle, xlabel="Number of Grid Points", ylabel="Error", label="", title="Error vs Number of Grid Points")
