
include("../../startup.jl")
using DifferentialEquations
using NumericallyIntegrateArrays
using Distributions
using LaTeXStrings
using LsqFit
using JLD2
using Interpolations

@load "radiotrack.jld2" t meansd semsd stdsd
t_rt = t
msd_rt = meansd;
sem_rt = semsd;
std_rt = stdsd;

function centreddiffusion2d!(du,u,p,t)
    D,h,n = p
    du[1] = (D/h^2)*(u[2]-u[1]) + (D/(2h^2))*(u[2]-u[1])
    for i = 2:n-1
        du[i] = (D/h^2)*(u[i-1]-2*u[i]+u[i+1]) + (D/(2i*h^2))*(u[i+1]-u[i-1])
    end
    du[n] = (D/h^2)*(u[n-1]-u[n]) + (D/(2n*h^2))*(u[n]-u[n-1])
end

function diffusion_sim2d(D, dt, u0, N, h, t, r)
    #us = [u0]
    tspan = (0.0, dt)
    msd = [2*pi*sum(u0.*r.^3)]
    us_int = [2*pi*sum(r.*u0)]
    us = [u0]
    for (i,ts) in enumerate(t[1:(end-1)])
        prob = ODEProblem(centreddiffusion2d!,u0,tspan,(D,h,N))
        u0 = solve(prob, Euler(),dt=dt).u[end]
        append!(msd,2*pi*sum(u0.*r.^3))
        append!(us,[u0])
        append!(us_int,2*pi*sum(u0.*r))
    end
    return msd, us_int, us
end

function shrink_rate(ts,R0 )
    a = R0^2 ./ts[end]^2
    Rt2 = R0^2 .- a.*ts.^2
    Rt2[Rt2.<0] .= 0.0
    return sqrt.(Rt2)
end

function domain_shrink2d(θ, D, dt,u0, N, h, t, tD, r)
    boxsize = zeros(Int, size(t))
    boxsize[t .<= θ[2]] .= div(θ[1],h)
    boxsize[t .> θ[2]] .= Int.(ceil.(shrink_rate(t[t .> θ[2]].-θ[2],θ[1])./h))
    boxsize[boxsize.>div(θ[1],h)] .= div(θ[1],h)
    boxsize[boxsize .== 1] .= 2
    boxsize[end] = 1
    tspan = (0.0, dt)
    msd = [2*pi*sum(u0.*r.^3)]
    us_int = [2*pi*sum(r.*u0)]
    for (i,ts) in enumerate(t[1:end-1])
        if ts > tD
            ns = length(u0)
            prob = ODEProblem(centreddiffusion2d!,u0,tspan,(D,h,ns))
            u0 = solve(prob, Euler(),dt=dt).u[end]
        end
        if ts > tD
            if boxsize[i+1] < length(u0) && boxsize[i+1] < boxsize[i]
                u0[boxsize[i+1]] += sum(u0[boxsize[i+1]+1:end])
                u0 = u0[1:boxsize[i+1]]
                r = r[1:boxsize[i+1]]
            end
        end
        append!(msd,2*pi*sum(u0.*r.^3))
        append!(us_int,2*pi*sum(u0.*r))
    end
    return msd, us_int
end

function r2_test(x1,y1,x_model,y_model)
    sstot = 0
    ssres = 0
    for (i,xi) in enumerate(x1)
        sstot += (y1[i]-mean(y1))^2
        ssres +=  (y_model[findfirst(x -> x==xi, x_model)]-y1[i])^2
    end
    return 1 - ssres/sstot
end


L = 2000
N = 50
h = L/N
D = 63.4#0.1h^2/2
u = zeros(N)
dt = (0.8h^2/2D)
T = t_rt[end]+dt
tD = 200.0

r = h:h:L
u0 = zeros(size(r))
u0[1] = 1/(2*pi*h)

t = dt:dt:T

sample_size = 10000

R_0s = rand(Uniform(1000,3000), sample_size)
t_ss = rand(Uniform(0.0,7.2e3), sample_size)


knots = (t,)
r2 = []

for i = 1:sample_size
    msds, u_int = domain_shrink2d([R_0s[i],t_ss[i]], D, dt,u0, N, h, t, tD, r);
    itp = interpolate(knots, msds, Gridded(Linear()))
    msd_itp = [ itp(ts) for ts in t_rt]
    append!(r2,r2_test(t_rt,msd_rt,t_rt,msd_itp))
end

@save "shrinkABC.jld2" R_0s, t_ss, r2
