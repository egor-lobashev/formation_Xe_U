include("rk4.jl")

using Printf

V = 3.6e5 # GLOBAL

function mid_size(C)
    return sum(Vector(1:120) .* C) / sum(C)
end

function print_results(trange,u)
    steps = Int((length(u) - 1) / 100)
    for i in [0,25,50,75,100] * steps
        if i == 0
            print("xj[", 1, "]=[")
        else
            print("xj[", Int(i/steps), "]=[")
        end
        for j in 1:120
            # @printf("%.0f, %10f\n", j, u[i][j]*V)
            @printf("%.6f", u[i+1][j]*V*20)
            if j != 120
                print(',')
            end
        end
        println(']')
    end

    # for i in 1:length(u)
    #     println(i, ", ", mid_size(u[i]))
    # end
end

function N_vac(n)
    return (2.65 - 2.24 * exp(-0.327 * n)) * n
end

function C_vac(C::Vector{Float64})
    ans = 0
    for n in 1:120
        ans += N_vac(n) * C[n]
    end
    return ans
end

C_sia_0 = 0.0002908 # GLOBAL
N = 120 # GLOBAL

function right(t::Float64, C::Vector{Float64})
    ans = Vector{Float64}(undef, N)
    ans .= 0.0

    C_sia = C_sia_0 + (C_vac(C) - C_vac_0)
    # C_sia = C_sia_0

    for i in 1:N
        for j in 1:i
            delta_C = factor[i, j] * C[i] * C[j] * C_sia
            # if delta_C*(V*20) > 0.5
            #     println(i, ", ", j, ", ", factor[i, j], ", ", C[i]*(V*20), ", ", C[j]*(V*20), ", ", C_sia*(V), "   =   ", delta_C*(V*20))
            # end

            ans[i] -= delta_C
            ans[j] -= delta_C

            if i+j <= N
                ans[i+j] += delta_C
            end
        end
    end

    return ans
end

include("matrix_D_sia.jl")
m11 = factor[1, 1]
# factor *= 2.5
# factor[1, 1] = m11

C_0 = Vector{Float64}(undef, N)
# for n in 1:length(C_0)
#     C_0[n] = 0.0005865 * exp(-1.3 * n)
# end
py_C_0 = [1157, 296, 111, 41, 26, 9, 9, 1, 2, 1, 1]
for n in 1:length(py_C_0)
    C_0[n] = py_C_0[n]/(V*20)
end
C_vac_0 = C_vac(C_0) # GLOBAL


problem = CauchyODEProblem(; f=right, tstart=0, tend=101, u₀=C_0)
trange, u = rk4(problem, nsteps=1000)

print_results(trange, u)