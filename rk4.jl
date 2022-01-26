struct CauchyODEProblem{T<:Real,F<:Function}
    bound::Tuple{T,T}   # отрезок интегрирования
    u₀::Vector{T}       # начальное значение интегрируемой функции
    f::F                # правая часть ОДУ
    function CauchyODEProblem(; f::Function, tstart::Real, tend::Real, u₀::Vector)
        new{Float64, typeof(f)}(
            float.(tuple(tstart, tend)),
            float.(u₀),
            f,
        )
    end
end

function rk4(problem::CauchyODEProblem; nsteps::Integer)
    u = Vector{Vector{Float64}}(undef, nsteps + 1)
    u[1] = problem.u₀
    tstart, tend = problem.bound
    trange = range(tstart, tend; length=nsteps+1)
    τ = step(trange)

    for i in 1:nsteps
        tᵢ, uᵢ = trange[i], u[i]
        
        k₁ = problem.f(tᵢ, uᵢ)
        k₂ = problem.f(tᵢ + τ/2, uᵢ + τ*k₁/2)
        k₃ = problem.f(tᵢ + τ/2, uᵢ + τ*k₂/2)
        k₄ = problem.f(tᵢ + τ, uᵢ + τ * k₃)

        u[i+1] = uᵢ + τ * (k₁ + 2*(k₂ + k₃) + k₄)/6
    end

    return trange, u
end