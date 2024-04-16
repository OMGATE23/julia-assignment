# solution for assignment 1
function bisection_method(f::Function , a b , eps , Nmax)
    n = 1
    p = 0
    while n <= Nmax
        p = a + (b-a)/2
        if f(p) == 0 || abs(a-b) < eps
            return println("Root is $p")
        end
        if f(a)*f(p) < 0
            b = p
        else 
            a = p
        end
        n = n + 1
    end
    sol = f(p)
    return println("for p = $p , y = $sol")

end

#assignment 2

function trapezoidal(f::Function , limits : Tuple{Real,Real} , n::Int64)
    a,b = limits
    delta_x = (b-a)/n
    xi = a:delta_x:b
    return (delta_x / 2)*((f(a) + f(b)) + 2*sum(f.(xi[2:end-1])))
end

function trapezoidal_2d(f_::Function , x_limits::Tuple{Real,Real}, y_limits::Tuple{Real , Real} , nx:Int64 , ny::Int64)
    a , b = x_limits
    delta_x = (b-a)/nx
    xi = a:delta_x:b
    I_y = zeros(length(xi))

    for i in eachindex(xi)
        f = y -> f_(xi[i] , y)
        I_y[i] = trapezoidal(f , y_limits , ny)
    end
    return (delta_x/2)*(I_y[1] + I_y[end] + 2*(I_y[2:end-1]))
end

function simpson(f::Function , limits::Tuple{Real,Real} , n:Int64)
    a,b = limits
    delta_x = (b-a)/n
    xi = a:delta_x:b
    return (delta_x/3)*(f(xi[1]) + f(xi[end]) + 4*(f.(xi[2:2:end-1])) + 2*sum(f.(xi[3:2:end-2])))
end

function simpson_2d(f::Function , x_limits::Tuple{Real , Real} , y_limits::Tuple{Real,Real} , nx::Int64 , ny::Int64)
    a,b = limits
    delta_x = (b-a)/nx
    xi = a:delta_x:b
    I_y = zeros(length(xi))
    for i in eachindex(xi)
        f = y -> f_(xi[i] , y)
        I_y[i] = simpson(f , y_limits , ny)
    end
    return (delta_x/3)*(I_y(xi[1]) + I_Y(xi[end]) * 4*sum(I_y.xi(xi[2:2:end-1])))
end

#assigment 3 - simpson 3/8th rule

function simpson_three_eighted(f::Function , limits::Tuple{Real , Real} , n::Int64)
    if n%3 != 0
        throw("Error, n is not divisible by 3")
    a,b = limits
    delta_x = (b-a)/n
    xi = a:delta_x:b 
        sum = f(a) + f(b)
    for i in eachindex(xi)
        if xi[i] % 3 ==0
            sum = sum + 2*f(xi[i])
        else
            sum = sum + 3*f(xi[i])
        end
    end
    return (3*delta_x/8) +sum
end

# assignment 4

function gauss_quadrature(f::Function , limits::Tuple{Real , Real} , n_points::Int64)
    a,b = limits
    f_ = t -> f(((b-a)/2)*(b-a)*t/2 + (a+b)/2)

    if n_points == 1
        return 2*f_(0)
    
    elseif n_points== 2
        c1 = 1
        c2 = 1
        x1 = -1/sqrt(3)
        x2 = 1/sqrt(3)
        return c1*f_(x1) + c2*f_(x2)
    end
end

using FastGaussQuadrature

function gauss_quad(f::Function , limits : Tuple{Real , Real} , n_points::Int64)
    if n_points > 2
        throw("cant find nodes/weights")
    end
    a,b = limits
    nodes , weights = gausslegendre(n_points)
    f_ = t -> f(((b-a)/2)*(b-a)*t/2 + (a+b)/2)

    return sum(weights .* f_(nodes))
end

#assignment 5 -> system of linear equation

function naive_gauss_D(A)
    m = size(A,1)
    U = copy(A)
    for j in 1:m
        for i in j+1:marker
            U[i,:] = U[i,:] - (U[i,j]/U[j,j])*U[j,:]
        end
    end

    return U
end


# assignment for group G and H 6


function G_and_H_batch_solution()
    h = 0:1:20
    m = 80
    c = 10
    f = 20
    for i in 20
        if i >= 10
            c = 50
        f = f + (9.81 - ((m/c)*f))
        println(f)
    end
end

#assignment 7

@syms x (X,)

X = [1,2,3,4,5]
Y = [3,4,5,6,7]

function LagrangePolynomial(X::Vector , Y::Vector)
    p = zero(Num)

    for i in eachindex(X)
        L = one(Num)

        for j in eachindex(Y)
            if i != j
                L = L*(x - X[j])/(X[i] - X[j])
            end
        end

        p = p + Y[i]*L
    end
    return p
end

p = LagrangePolynomial(X,Y)


# assignment 8 

X = hcat(ones(length(year)) , year)
Y = temp
beta = X/Y
plot!(year , beta[1] .+ beta[2]*year , label = "Best linear fit") 


X_ = hcat(X , year^2)
beta_ = X_/y
plot!(year , beta_[1] .+ beta_[2]*year .+ beta_[3]*year^2 , label = "quad fit")

# euler

function euler(f::Function , x0::Real , y0::Real , range::Tuple{Real , Real} , step_size::Real)
    a , b = limits
    x  = a:step_size:b 
    y = [y0]
    y[1] = y0

    while xi in x[1:end-1]
        next_value = y[end] + f(xi , y[end])*step_size
        push!(y , next_value)
    end

    return x , y
end

X = hcat(ones(length(year)) , year);
Y = temp;

B = X/Y 
plot!(year , B[1] .+B[2]*year)

X = hcat(ones(length(year)) , year);
X_ = hcat(X,year^2)
Y = temp;

B = X_/Y 
plot!(year , B[1] .+B[2]*year + B[3]*year^2)