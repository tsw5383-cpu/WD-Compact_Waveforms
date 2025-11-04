using Documenter, Integrate

makedocs(modules = [Integrate], sitename = "WD_NS_Interaction.jl",
         format = Documenter.HTML(), checkdocs = :exports)

deploydocs(repo = "github.com/PsuAstro528/project-wd_ns_interaction.git", target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true,
)
    # , devbranch="main")
