using LinearAlgebra



A = [0.02 0.01 0 0; 1 2.0 1 0.0; 0 1 2 1; 0 0 100 200]


b = [0.02; 1; 4; 800]

A_ = [A b]

function gausselimiation_with_partialpivot(A_::Matrix)
    m = size(A_, 1)
    U = copy(A_)
    for j in 1:m             
        max_val, max_ind =  findmax(abs.(U[j:end, j]))
        if max_val ≠ U[j,j]
            U[j,:] , U[max_ind + (j-1),:] = U[max_ind + (j-1),:], U[j,:]
        end

        for i in j+1:m
            if U[i,j] == 0.0
                continue
            else
                U[i,:] = U[i,:] - (U[i,j]/U[j,j])*U[j,:]
            end
        end
    end

    return U
end

function backsubstitution(A::Matrix,b::Vector)
    m = size(A, 1)
    x = similar(b)

    for j= m:-1:1
        r = b[j]
        for k = j+1:m  # this indicates all the previous elements
            r = r  - (A[j,k]*x[k])  # x[k]  is previously found solution

        end
        x[j] = r/A[j,j]
    end
    return x
end
U = gausselimiation_with_partialpivot(A_)
x = backsubstitution(U[:,1:end-1], U[:,end])

println("solution" , x)



# using LinearAlgebra
# """
#    .... + A[i,i-1]*x[i-1] + A [i,i]*x[i] + A[i, i+1]*x[i+1] +... = b[i]

#    x[i] = (b[i] - (...+ A[i,i-1]*x[i-1] + A[i, i+1]*x[i+1]  +.... ) )/A [i,i]

#    dot product of  A[i,all indices except i], x[all indices except i]

# """

# function jacobi(A::Matrix, b::Vector, x_old::Vector, tol::Float64)
#     @assert size(b) ==  size(x0)

#     err = norm(A*x_old - b)  # means         || Ax -b ||

#     x_next = similar(x_old)
#     while err > tol
#         println("error is ",err)
#         for i in eachindex(b)  #same  as 1: length(b)... i is for each equation
#             sum_ = dot(A[i, 1:end .≠ i ], x_old[1:end .≠ i])   # look here

#             x_next[i] = (b[i] - sum_ )/A[i,i]
#         end
#         x_old = x_next  # updating next guess as old in jacobi
#         err = norm(A*x_old - b)
#     end

#     return x_next
# end

# # new example...

# A = rand(50,50) + 10*I(50)
# b = rand(50)
# x0 = zeros(50)
# tol=1e-4
# sol  = jacobi(A,b, x0, tol)



# function gaussSiedal(A::Matrix, b::Vector, x_old::Vector, tol::Float64)
#     @assert size(b) ==  size(x0)

#     err = norm(A*x_old - b)  # means         || Ax -b |

#    # x_next = similar(x_old)

#     while err > tol
#         for i in eachindex(b)
#             sum_ = dot(A[i, 1:end .≠ i ], x_old[1:end .≠ i])   # look here
#             x_old[i] = (b[i] - sum_ )/A[i,i]
#         end
#         err = norm(A*x_old - b)  
#     end

#     return x_old
# end
# sol  = gaussSiedal(A, b, x0, tol)

# function gauss_quad(f::Function, limits::Tuple{Real , Real} , n_points::Int64)
#     if n_points > 2
#         throw("cant find weights for n>2")
#     end
#     a,b = limits

#     f_ = t -> ((b-a)/2)*f((b-a)*t/2 + (a+b)/2)

#     if n_points == 1
#         println("answer")
#         return 2*f_(0)

#     elseif n_points == 2
#         c1 = 1;
#         c2 = 1;
#         x1 = -1/sqrt(3);
#         x2 = 1/sqrt(3)
#         println("answer")
#         return c1*f_(x1) + c2*f_(x2)
#     end
    
# end

# using FastGaussQuadrature


# function gauss_quad_smart(f::Function, limits::Tuple{Real , Real} , n_points::Int64)
#     if n_points > 2
#         throw("cant find weights for n>2")
#     end
#     a,b = limits
#     nodes, weights = guasslegendre(n_points)
#     f_ = t -> ((b-a)/2)*f((b-a)*t/2 + (a+b)/2)

#     return sum(weights .* f_(nodes)) 
    
# end

# f = x -> ((sin(x))^2)

# gauss_quad(f , (0.0 , pi) , 2)
# gauss_quad_smart(f , (0.0 , pi) , 2)

# function trapezoidal(f::Function , limits::Tuple{Real , Real} , n::Int64)
#     a,b = limits
#     delta_x = (b-a)/n_points
#     xi = a:delta_x:b
#     return (delta_x/2)*(f(xi[1]) + f(xi[end] + 2 * sum(f.(xi[2:end]))))
# end


