module StochasticSimulations
export diffusion_simulation, square_distance

using Distributions, Random, StatsBase

function diffusion_step(args,x,y)
    """
    Takes position x,y and returns a new position after a diffusion step
    with diffusion coefficient D and timestep dt
    """
    return x .+ rand(Normal(0,sqrt(2*args.D*args.dt)),args.N), y .+ rand(Normal(0,sqrt(2*args.D*args.dt)),args.N)
end

function diffusion_simulation(args)
    """
    Stochastic diffusion simulation with discrete time and arguments
    args = (D=diffusion coefficient,N=number of particles,T=time of simulation,dt=timestep)
    returns positions in the trajectory.
    """
    t = args.dt:args.dt:args.T
    x = zeros(Float64,args.N,length(t))
    y = zeros(Float64,args.N,length(t))
    for (j,ts) in enumerate(t[1:end-1])
        x[:,j+1],y[:,j+1]=diffusion_step(args,x[:,j],y[:,j])
    end
    return (x=x,y=x)
end

function square_distance(xy)
    """
    Calculates the squared distance from 0 for trajectory xy.
    Returns a tuple of mean squared distance and standard error.
    """
    SD = xy.x.^2 .+ xy.y.^2;
    SEM = Array{Float64}(undef,size(SD)[2])
    for i = 1:size(SD)[2]
        SEM[i] = sem(SD[:,i])
    end
    return (MSD=mean(SD,dims=1),std_error = SEM)
end


end
