module GW
    using NFFT, FFTW

    const Mpc_to_km_factor = 2*G*3.2e-20/c^4 #for conversion from Mpc to km in quadrupole formula
    
    """Given a luminosity distance `DL` in units of Mpc and an orientation parameterized by an inclination `i` and azimuthal angle `φ`, returns a function which calculates the instantaneous gravitational strain waveform based on the quadrupole formula."""
    function gravitational_waveform_ftn(inc::Float64, φ::Float64, DL::Float64; mode::String="Debug")::Function
        #fixes a function to return gravitational wave given orientation and distance
        #DL is in Mpc
        if mode=="Debug"
            @assert 0<=inc<=π
            @assert 0<=φ<2π
            @assert DL<10000 #Ensure correct units for DL
        end
        cosi = cos(inc)
        cosφ = cos(φ)
        sini = sin(inc)
        sinφ = sin(φ)
        e1=[cosi*cosφ, sini*sinφ]
        e2=[-sinφ,  cosφ]
        es=[e1,e2]
    
        matrices=zeros(4,2,2)
        for i in 1:2
            for j in 1:2
                @inbounds matrices[2*(i-1)+j ,:,:] = es[i] * es[j]' #OUTER PRODUCT OF VECTORS
            end
        end
    
        #given fixed outer products, calculate gravitational wave strain tensor
        function ftn(ddI::Matrix{Float64})::Matrix{Float64}
            prod=zeros(2,2) #will end up being in transverse gauge
            for i in 1:2
                for j in 1:2
                    @inbounds prod[i,j] = sum(matrices[2*(i-1)+j ,:,:].*ddI) #OUTER PRODUCT OF VECTORS    
                end
            end
            prod_TT = prod .- (prod[1,1]+prod[2,2]) .* [1/2 0; 0 1/2] #traceless
            return((Mpc_to_km_factor/DL).*prod_TT)
        end
    
        return ftn 
    end

    """ Calculate PSD for an unevenly spaced time series """
    function PSD(ts::Vector{Float64},y::Vector{Float64})
        @assert length(ts)==length(y)
    
        #normalizes times to [-0.5,+0.5] and calculate frequencies
        t_norm=(ts-ts[1])/(ts[end]-ts[1]) .- 0.5
        
        #instantiate and solve NFFT problem
        N = 2^ceil(Int, log2(length(ts)))
        #number of frequency bins to use
        problem = nfft_plan(t_norm, N)
        problem.y=y
        nfft_execute!(problem)

        #Calculate and return PSD
        problem.y .* conj(problem.y) / (ts[end]-ts[1])
    end
    

    """Given a gravitational waveform, calculates the characteristic strain of the plus and cross polarizations"""
    function characteristic_strain(waveform::Matrix{Float64})
        @assert length(waveform[1])==3
    
        ts=waveform[:,1] #time series
        hp=waveform[:,2] #plus polarization strain
        hc=waveform[:,3] #times polarization strain

        #PSDs for plus and cross polarizations
        PSD_p=PSD(ts,hp)
        PSD_c=PSD(ts,hc)

        #calculate output frequencies
        N = 2^ceil(Int, log2(length(ts)))
        freqs=FFT.fftfreq(N,mean(ts[2:end].-ts[1:end-1]))

        #characteristic strain
        hc_p = sqrt.(2 .* PSD_p)
        hc_c = sqrt.(2 .* PSD_c)

        #return output values
        out=zeros(N,3)

        out[:,1].=freqs
        out[:,2].=sqrt.(2.*PSD_p) #conventional definition for characteristic strain
        out[:,3].=sqrt.(2.*PSD_c)

        out
    end
    
end