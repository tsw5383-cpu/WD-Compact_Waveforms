using Pkg
Pkg.activate(".")
using BenchmarkTools
using Integrate

sim_conditions = (M_WDs= 0.5:0.01:0.8, M_NSs= 1.2:0.1:1.4)
## UNCOMMENT below lines of code and then add `optional_conditions=optional_conditions` at the end of the many_integrate function call parameter arguments to test running with optional_conditions.
## Example: Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,uneven=false,optional_conditions=optional_conditions)
#optional_conditions = (mu_es=(2.0,2.35), WD_mass_radius_equations=(1,2), As=(3,10,30),
#    use_polytropic_density=(true,false), a_ring_fracs=(0.9,1.0), 
#    mass_loss_modes=('j','i','c'), eta_accs=(0.1,0.2), use_period_calc=(true,false))

step_limit::Int=10000
#println("Inputted M_WDs: $(sim_conditions.M_WDs). Inputted M_NSs: $(sim_conditions.M_NSs).")
#distributed_parallel_time_many = @belapsed Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,uneven=false)
#println("The time for many_integrate() in the shared memory parallel version is $(distributed_parallel_time_many) sec.")

M_WD = 1.0
Integrate.integrate(M_WD, step_limit=step_limit,write_file=true,uneven=false)
#parallel_time = @belapsed Integrate.integrate($M_WD, step_limit=$step_limit,write_file=false,uneven=false)
#println("The time for integrate() in the shared memory parallel version is $(parallel_time) sec.")