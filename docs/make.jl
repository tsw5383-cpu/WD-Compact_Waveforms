using Documenter, Math, Integrate, Save

makedocs(modules = [Math, Integrate, Save], sitename = "WD_NS_Interaction.jl",
         format = Documenter.HTML())

deploydocs(repo = "github.com/PsuAstro528/project-template.git", target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true,
)
    # , devbranch="main")
