module Save
    using CSV, Dates
    
    """Formats the filename as YYYYMMDD_M_WD=*M_WD*_M_NS=*M_NS*_A=*A*.csv
		if a filename is not provided"""
	function format_filename(M_WD::Float64,M_NS::Float64, A::Float64;filename::Union{String,Int}=-1)::String
		#
        M_WD=round(M_WD,3)
        M_NS=round(M_NS,3)
        A=round(A,3)
		if filename==-1 
		        today_str = Dates.format(Dates.now(), "yyyymmdd")
		        base_filename = "Outputs_2/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_A=$(A).csv"
		        filename = base_filename
		        suffix = 2
		        while isfile(filename)
		            filename = "Outputs_2/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_A=$(A)_$(suffix).csv"
		            suffix += 1
				end
		end
		return filename
	end
		
	"""Writes array of format given in Integrate.integrate() to a file."""
	function save_as_csv(datumses::AbstractMatrix{<:Real},
					  M_WD::Float64,M_NS::Float64,A::Float64; filename::Union{String,Int}=-1)::Nothing
		
		
		filename=format_filename(M_WD,M_NS;filename=filename)
		# Writes the file to a CSV titled filename

		isdir("Outputs_2") || mkdir("Outputs_2")
		#creates output directory if not already extant
		
			CSV.write(filename, 
						(t = datumses[:,1],
     					ddI_xx = datumses[:,2],
     					ddI_xy = datumses[:,3],
						a = datumses[:,4],
     					M_WD = datumses[:,5],
     					M_NS = datumses[:,6],
						R_WD = datumses[:,7],
     					J = datumses[:,8],
     					θ = datumses[:,9]))
		return nothing
	end
end