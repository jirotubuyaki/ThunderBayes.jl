using Documenter
using ThunderBayes

makedocs(
    sitename = "ThunderBayes",
    format = Documenter.HTML(),
    modules = [ThunderBayes]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
